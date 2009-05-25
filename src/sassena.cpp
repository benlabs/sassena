/*
 *  sassena.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "common.hpp"

// standard header
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <libconfig.h++>
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "sample.hpp"
#include "scatterdata.hpp"
#include "settings.hpp"
#include "tasks.hpp"

using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace boost::filesystem;
namespace mpi = boost::mpi;

// user-defined MPI TAGS
enum MPI_TAGS {
	MPI_TAG_NOTHING,
	MPI_TAG_READY,
	MPI_TAG_HANGUP,
	MPI_TAG_CLOSING,
	MPI_TAG_MSG,
	MPI_TAG_IDLE,	
	
	MPI_TAG_SCATTER,
	MPI_TAG_SCATTERRESULT
};

enum MPI_VALUES {
	MPI_HANGUP
};

//forward
void main_slave(boost::mpi::environment& env, boost::mpi::communicator& world,int argc,char** argv);

void print_eta(int percent,int progress,timeval start,timeval end) {
		gettimeofday(&end, 0);
		int seconds = (end.tv_sec-start.tv_sec);
		int microseconds = (end.tv_usec-start.tv_usec);
		microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
		seconds = (microseconds < 0) ? seconds-1 : seconds;
		double etotal = 100 * ( seconds + (microseconds/1000000.0) ) / percent;
		double eta = etotal - seconds + (microseconds/1000000.0);
		clog << "INFO>> " << "Progress: " << percent << "%" << " , ETOTAL(s): " <<  etotal  << " , ETA(s): " <<  eta  << endl;
}

vector<int> get_step_assignments(size_t thisrank,size_t nodes, size_t frames) {
	
	vector<int> result; 
	size_t current_step_size=0;
	size_t current_node = nodes-1;
	int step_var=1; // go forward
	while(current_step_size<=frames) {
		if (thisrank==current_node) result.push_back(current_step_size);
		current_node += step_var;					
		if ( (current_node<0) ) {
			current_node=0;
			step_var *= -1.0; // reflect
		} 
		if (current_node>nodes-1) {
			current_node=nodes-1;
			step_var *= -1.0; // reflect
		}
		current_step_size++;
	}
	
	return result;
}

int main(int argc,char** argv) {

	//------------------------------------------//
	//
	// MPI Initialization
	//
	//------------------------------------------//	
	
  	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;
    // Get the number of processors this job is using:
	int rank = world.rank();

    // Get the rank of the processor this thread is running on.  (Each
    // processor has a unique rank.)
	int size = world.size();
      
    // The rank 0 node is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:
    // cerr << "ERROR>> (" << rank << "/" << size ")" << " some error message " << endl;
	// clog << "INFO>> " << "############ This is a debug hello from node ############             " << rank << endl;

	Sample sample;

	if (rank==0) {
		//------------------------------------------//
		//
		// Some welcome message....
		//
		//------------------------------------------//	
	
		clog << "INFO>> " << "This software is being developed by Benjamin Lindner and Franci Merzel.  " << endl;
		clog << "INFO>> " << "For help, suggestions or correspondense use:                             " << endl;
		clog << "INFO>> " << "ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  " << endl;
		clog << "INFO>> " << "franc@cmm.ki.si, Franci Merzel (Main Developer, Methodology )            " << endl;
		clog << "INFO>> " << "For publications use the following references:                           " << endl;
		clog << "INFO>> " << "........................................................................." << endl;
		clog << "INFO>> " << "1. SASSIM: a method for calculating small-angle X-ray and                " << endl;
		clog << "INFO>> " << "   neutron scattering and the associated molecular envelope from         " << endl;
		clog << "INFO>> " << "   explicit-atom models of solvated proteins, F. Merzel and J. C. Smith, " << endl;
		clog << "INFO>> " << "   Acta Cryst. (2002). D58, 242-249                                      " << endl;	
		clog << "INFO>> " << "........................................................................." << endl;
		clog << "INFO>> " << "SASSIM is a software for calculating scattering spectra from all-atomic.." << endl;
		clog << "INFO>> " << "" << endl;

		//------------------------------------------//
		//
		// Test the mpi environment
		//
		//------------------------------------------//	
	
		//------------------------------------------//
		//
		// Setting up the sample
		//
		//------------------------------------------//	
	
		clog << "INFO>> " << "Reading configuration file " << argv[1] << endl;
		// read in configuration file, delete command line arguments
		if (!Settings::read(argc,argv)) { cerr << "Error reading the configuration" << endl; throw; }

	   clog << "INFO>> " << "Reading structure from: " << (const char *) Settings::get("main")["sample"]["structure"]["file"] << endl;

	   // create the sample via structure file	
		sample.add_atoms(Settings::get_filepath(Settings::get("main")["sample"]["structure"]["file"]), Settings::get("main")["sample"]["structure"]["format"]);

	   // add selections / groups
	   libconfig::Setting* sgroups = &Settings::get("main")["sample"]["groups"];
	   for (int i=0;i<sgroups->getLength();i++) {

			string gn = (*sgroups)[i].getName(); // groupname
			string fn = (*sgroups)[i]["file"]; // filename
			string ft = (*sgroups)[i]["type"]; // filetype
			if ((*sgroups)[i].exists("select")) {
				string sn = (*sgroups)[i]["select"]; // selection field name
				double sv = (*sgroups)[i]["select_value"]; // selection field value, positive match			
				sample.add_selection(gn,fn,ft,sn,sv);
			}
			else {
				sample.add_selection(gn,fn,ft);			
			}

			clog << "INFO>> " << "Adding selection: " << gn << endl;	
		}
	
		// add custom selections here ( if not already set! )
		if (sample.atomselections.find("system")==sample.atomselections.end()) {
			// this shortcut creates a full selection
			sample.add_selection("system",true);		
		}
	
		// apply deuteration
		if (Settings::get("main")["sample"].exists("deuter")) {
			if (Settings::get("main")["sample"]["deuter"].getType()==libconfig::Setting::TypeString) {
				clog << "INFO>> " << "Deuteration of group " << (const char *) Settings::get("main")["sample"]["deuter"] << endl;		
				sample.deuter(Settings::get("main")["sample"]["deuter"]);
			}
			else if (Settings::get("main")["sample"]["deuter"].getType()==libconfig::Setting::TypeList) {
		    	for (int i=0;i<Settings::get("main")["sample"]["deuter"].getLength();i++) {
					clog << "INFO>> " << "Deuteration of group " << (const char *) Settings::get("main")["sample"]["deuter"][i] << endl;		
					sample.deuter(Settings::get("main")["sample"]["deuter"][i]);
				}
			}
		}
		else {
			clog << "INFO>> " << "No explicit deuteration" << endl;		
		}
	
		// read in frame information
		for (int i=0;i<Settings::get("main")["sample"]["frames"].getLength();i++) {
			clog << "INFO>> " << "Reading frames from: " << (const char *) Settings::get("main")["sample"]["frames"][i]["file"] << endl;
	        size_t nof = sample.frames.add_frameset(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"],sample.atoms);
			clog << "INFO>> " << "Found " << nof << " frames" << endl;			
	//		ifstream dcdfilelist(Settings::get_filepath(dcdfilename).c_str());
	//		std::string line;
	//
	//		while (getline(dcdfilelist,line)) {
	//			add_file(line,atoms,(recursion_trigger-1));
	//		}
	//		sample.add_frame(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"]);
		}
		// final statement: total number of frames
		clog << "INFO>> " << "Total number of frames found: " << sample.frames.size() << endl;			

		// select wrapping behaviour
		if (Settings::get("main")["sample"]["pbc"]["wrapping"]) {
			sample.frames.wrapping = true;
			string centergroup = (const char *) Settings::get("main")["sample"]["pbc"]["center"];
			sample.frames.centergroup_selection = sample.atomselections[centergroup];
			clog << "INFO>> " << "Turned wrapping ON with center group " << centergroup << endl;				
		}
		else {
			sample.frames.wrapping = false;		
			clog << "INFO>> " << "Turned wrapping OFF "  << endl;						
		}



		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		string m = Settings::get("main")["scattering"]["background"]["method"];
		bool rescale = Settings::get("main")["scattering"]["background"]["rescale"];
		// rescale automatically triggers background analysis routine. Maybe independent in the future...
		if (m=="auto" || ((m == "fixed") && rescale)) {
			if (m=="fixed") clog << "INFO>> " << " rescale=true triggered analysis of the background scattering factor!" << endl;

			// determine background scattering density from grid based analysis
			// necessary:
			// group = name of the group which accounts for the background
			// resolution = grid resolution
			// hydration = hydration layer around each particle, 
			//    ie. solvent molecules within this grid distance are not counted towards bulk
			int r     = Settings::get("main")["scattering"]["background"]["resolution"]; // resolution of grid
			double h  = Settings::get("main")["scattering"]["background"]["hydration"]; // hydration layer of each particle (not background)
		
			clog << "INFO>> " << "Calculate background scattering length density" << endl;
			clog << "INFO>> " << "using grid resolution " << r << " and hydration of " << h << endl;		
			pair<double,double>	b0 = Analysis::background_avg(sample,r,h,CartesianCoor3D(0,0,0));
			sample.background = b0.first;
			clog << "INFO>> " << "Using average background scattering length: " << b0.first << " +- " << b0.second << endl;		
		}
		else if (m=="fixed") {
			double f = Settings::get("main")["scattering"]["background"]["value"];
			sample.background = f;
			clog << "INFO>> " << "Setting background scattering length density to fixed value of " << f << endl;
		}
		else if (m=="none") {
			// switch off background correction, i.e. assume system in vaccuum
			sample.background = 0.0;
		}

		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

		sample.frames.clear_cache(); // reduce overhead

		clog << "INFO>> " << "Exchanging sample information with compute nodes... " << endl;

		world.barrier();

		broadcast(world,sample,0);
	}
	else {

		if (!Settings::read(argc,argv)) { 
			cerr << "ERROR>> (" << rank << "/" << size << ")" << " Error reading the configuration" << endl;
			return -1;
		}

		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, sample);	
	}
	

	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//

	// for timing:
 	timeval start, end;	

	// prepare an array for all qqqvectors here:
	std::vector<CartesianCoor3D> qqqvectors;
	std::vector<int> frames;
	
	// fill qqqvectors
	Settings::get_qqqvectors(qqqvectors);
	
	// fill frames
	int framestride = 1;
	if (Settings::get("main")["scattering"].exists("framestride")) {
		framestride = Settings::get("main")["scattering"]["framestride"];
	}
	for(size_t i=0;i<sample.frames.size();i+= framestride) {
		frames.push_back(i);
	}	

	int taskcounter = 0; int progress = 0;
	
	// generate Tasks
	string mode;
	if (Settings::get("main")["scattering"].exists("correlation")) {
		mode = (const char * ) Settings::get("main")["scattering"]["correlation"]["type"];
	}
	else {
		mode = "none";
	}

	Tasks tasks(frames,qqqvectors,world.size(),mode);
	if (rank==0) { 
		clog << "INFO>> " << "Task assignment:" << endl;		
		tasks.print(); 
	}
	
	Tasks my_tasks = tasks.scope(rank);
		
	// get color(current worker class ID) from first task
	// fix this for cases where my_tasks is empty!
	boost::mpi::communicator local = world.split(my_tasks[0].color);
	
//	map<CartesianCoor3D<map<size_t,double> > keyvals;
	string target = Settings::get("main")["scattering"]["target"];

	clog << "INFO>> " << "Scattering target selection: " << target << endl;

	gettimeofday(&start, 0);
	int progess = 0;
	
	map<pair<CartesianCoor3D,size_t>, double > qF_Task_results; // used if mode = "none"
	map<pair<CartesianCoor3D,size_t>, double > q_Task_results; // used if mode = "time"
	
	for (Tasks::iterator ti=my_tasks.begin();ti!=my_tasks.end();ti++) {

			// progress indicator
			if (rank==0) {
				int percent = int(taskcounter*100/my_tasks.size())  ;
				if ( percent >= progress  ) {
					print_eta(percent,progess,start,end);
					if (percent<10) {
						progress = 1*(percent/1) + 1;									
					}
					else {
						progress = 10*(percent/10) + 10;									
					}
				}				
			}

			// set scattering amplitudes for current q vector
			Analysis::set_scatteramp(sample,sample.atomselections[target],ti->q,true);
			
			// generate q-vector list for orientational averaging
			vector<CartesianCoor3D> qvectors;

			if (Settings::get("main")["scattering"].exists("average")) {
			
				libconfig::Setting& s = Settings::get("main")["scattering"]["average"];
				double resolution = -1.0;
				if (s.exists("resolution")) {
					resolution = s["resolution"];
				} 

				string avtype = s["type"]; 
			
				uint32_t qseed = 1; // ti->qseed			
				if (avtype=="none") {
					qvectors.push_back(ti->q);
				}
				else if (avtype=="sphere") {
					string avmethod = s["method"];
					if (avmethod=="bruteforce") {
						if (resolution==-1.0) resolution=1.0;

						string avvectors = s["vectors"];
						Analysis::qvectors_unfold_sphere(avvectors,ti->q,qseed,resolution,qvectors);
					}
					if (avmethod=="multipole") {
						// don't unfold. 
						// multipole expansion works on one qvector
					}		
				}
				else if (avtype=="cylinder") {
					string avmethod = s["method"];		
					if (avmethod=="bruteforce") {
						if (resolution==-1.0) resolution=1.0;
						string avvectors = s["vectors"];
						Analysis::qvectors_unfold_cylinder(avvectors,ti->q,qseed,resolution,qvectors);
					}
					if (avmethod=="multipole") {
						// don't unfold. 
						// multipole expansion works on one qvector
					}
				}	
				else {
					cerr << "ERROR>> " << "unrecognized averaging type. Use 'none' if no averaging is to be done" << endl;
					throw;
				}
			
			}
			else {
				qvectors.push_back(ti->q);
			}

			// block qqqvectors;
			size_t qvector_blocking = 10;
			
//			for(size_t qvector_block = 0; qvector_block <= (qvectors.size()/qvector_blocking); ++qvector_block) {
//				size_t qii = 0;
//				vector<CartesianCoor3D> qvectors_sub;
//				while ( (qvector_blocking*qvector_block + qii)<qvectors.size() && (qii<qvector_blocking) ) {
//					qvectors_sub.push_back( qvectors.at(qvector_blocking*qvector_block + qii) );
//					qii++;
//				}
//				// first element: qvectors[qvector_blocking*qvector_block + i ]

				vector<int> frames_i = ti->frames(rank);
				std::string interference_type = Settings::get("main")["scattering"]["interference"]["type"];
				if (interference_type == "self") {

					map<int,vector<vector<complex<double> > > > scatbyframe; // frame -> qvectors/atoms/amplitude
			
					// iterate through all frames this node is supposed to do
					for(size_t i = 0; i < frames_i.size(); ++i)
					{
						Atomselection& as = sample.atomselections[target];
						sample.frames.load(frames_i[i],sample.atoms,sample.atomselections[target]);				

						// holds the scattering amplitudes for the current frame:
						vector<vector<complex<double> > >& scattering_amplitudes = scatbyframe[frames_i[i]];

						libconfig::Setting& s = Settings::get("main")["scattering"]["average"];
						string avtype = s["type"]; 						
						if (avtype=="none") {
							// qseed not used here
							vector<complex<double> > scat;
							Analysis::scatter_none(sample,as,ti->q,scat);
							scattering_amplitudes.push_back(scat);
						}
						else if (avtype=="sphere") {
							string avmethod = s["method"];
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,as,qvectors,scattering_amplitudes);						
							}
							if (avmethod=="multipole") {
//								if (resolution==-1.0) resolution=17.0;
								throw;
//					   			Analysis::scatter_sphere_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
							}		
						}
						else if (avtype=="cylinder") {
							string avmethod = s["method"];		
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,as,qvectors,scattering_amplitudes);		
							}
							if (avmethod=="multipole") {
//								if (resolution==-1.0) resolution=10.0;
								throw;
//					   			Analysis::scatter_cylinder_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
							}
						}	
						else {
							cerr << "ERROR>> " << "unrecognized averaging type. Use 'none' if no averaging is to be done" << endl;
							throw;
						}	
					} // frames processed

					// aggregate results
					if (ti->mode=="none") {

						// no correlation -> instanenous scattering intensity
						for(map<int,vector<vector<complex<double> > > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++) {
							double scatsum=0;				
							for(vector<vector<complex<double> > >::iterator si=sbfi->second.begin();si!=sbfi->second.end();si++) {
								for(vector<complex<double> >::iterator si2=si->begin();si2!=si->end();si2++) {							
									scatsum += abs(conj(*si2)*(*si2));
								}
							}
							qF_Task_results[make_pair(ti->q,sbfi->first)] = scatsum/(sbfi->second.size()); 
						}
					}
					else if (mode=="time") {

						// check size of vector<scattering amplitudes>, this is the size of unfolded q vectors
						size_t aqscount = 0;
						for(map<int,vector<vector<complex<double> > > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++) {
							if (aqscount==0) aqscount = sbfi->second.size(); else if (aqscount!=sbfi->second.size()) throw;
						}

						for(size_t asqi = 0; asqi < aqscount; ++asqi)
						{
						// communicate vector of current aqs

							// sbfi->first = frame; sbfi->second = qvector/atoms/scattering-amplitudes
							vector<pair<int,vector<complex<double> > > > aq_vectors_by_frame;				
							for(map<int,vector<vector<complex<double> > > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++)
							{
								aq_vectors_by_frame.push_back(make_pair(sbfi->first,sbfi->second[asqi]));
							}
							vector<vector<pair<int,vector<complex<double> > > > > vector_in,vector_out;
							vector_in.resize(local.size(),aq_vectors_by_frame);

							boost::mpi::all_to_all(local,vector_in,vector_out);

							// decompose vector_out into new table: frame <-> Aqs, use frame number as implicit position
							vector<vector<complex<double> > > aq_vectors; aq_vectors.resize(frames.size());
							for(size_t i = 0; i < vector_out.size(); ++i)
							{
								for(size_t j = 0; j < vector_out[i].size(); ++j)
								{
									aq_vectors[vector_out[i][j].first]=vector_out[i][j].second;
								}
							}

							// now determine which correlation 'step' WE need to do:
							vector<int> stepsizes = get_step_assignments(local.rank(),local.size(),frames.size()/2); // frames.size()/2 is the length of the correlation we are interested in...

							map<int,complex<double> > my_AAconj;

							for(vector<int>::iterator ssi=stepsizes.begin();ssi!=stepsizes.end();ssi++) {
								complex<double> AAconj_sum=0;
								int AAconj_count=0;					
								size_t current_frame=0;
								while( (current_frame+(*ssi))<frames.size()) {
									// do a vector-vector multiply here:
									for(size_t i = 0; i < aq_vectors[current_frame].size(); ++i)
									{
										AAconj_sum += aq_vectors[current_frame][i]*conj(aq_vectors[current_frame+(*ssi)][i]);
									}
									// vector-vector multiply finished.
									AAconj_count++;
									current_frame++;
								}
								my_AAconj[*ssi] = AAconj_sum / double(AAconj_count);
							}

							vector<map<int,complex<double> > > all_AAconj;
							boost::mpi::gather(local,my_AAconj,all_AAconj,0);

							if (local.rank()==0) {
								// reorder results, use stepsize as implicit index
			//					vector<complex<double> > AAconj_series; AAconj_series.resize((frames.size()/2)+((frames.size()+1)%2)); // correlation length = size

								vector<complex<double> > AAconj_series; AAconj_series.resize((frames.size()/2)+1); // correlation length = size					
								for(vector<map<int,complex<double> > >::iterator ai=all_AAconj.begin();ai!=all_AAconj.end();ai++) {
									for(map<int,complex<double> >::iterator aii=ai->begin();aii!=ai->end();aii++) {
										AAconj_series[aii->first]=aii->second;
									}
								}
			//					ofstream ofile("test.dat",ios_base::app);
								for(size_t i = 0; i < AAconj_series.size(); ++i)
								{
									q_Task_results[make_pair(ti->q,i)] = AAconj_series[i].real(); 
									// check THIS! take the real part or do a conj-multiply?
			//						ofile << ti->q.x << "\t" << ti->q.y << "t" << ti->q.z << "\t" << i << "\t" << AAconj_series[i].real() << endl;
								}
			//					ofile << endl;
							}

						}

						// wait for all nodes to finish, before communicating
						local.barrier();
					}

					
				}
				else if (interference_type == "all") {

					map<int,vector<complex<double> > > scatbyframe; // frame -> qvectors/amplitude
			
					// iterate through all frames this node is supposed to do
					for(size_t i = 0; i < frames_i.size(); ++i)
					{
						sample.frames.load(frames_i[i],sample.atoms,sample.atomselections[target]);				

						Atomselection& as = sample.atomselections[target];
						// holds the scattering amplitudes for the current frame:
						vector<complex<double> >& scattering_amplitudes = scatbyframe[frames_i[i]];
						
						libconfig::Setting& s = Settings::get("main")["scattering"]["average"];
						double resolution = -1.0;
						if (s.exists("resolution")) {
							resolution = s["resolution"];
						}						
						string avtype = s["type"]; 							
						if (avtype=="none") {
							// qseed not used here
							complex<double> scat = Analysis::scatter_none(sample,as,ti->q);
							scattering_amplitudes.push_back(scat);
						}
						else if (avtype=="sphere") {
							string avmethod = s["method"];
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,as,qvectors,scattering_amplitudes);						
							}
							if (avmethod=="multipole") {
								if (resolution==-1.0) resolution=17.0;
					   			Analysis::scatter_sphere_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
							}		
						}
						else if (avtype=="cylinder") {
							string avmethod = s["method"];		
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,as,qvectors,scattering_amplitudes);		
							}
							if (avmethod=="multipole") {
								if (resolution==-1.0) resolution=10.0;
					   			Analysis::scatter_cylinder_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
							}
						}	
						else {
							cerr << "ERROR>> " << "unrecognized averaging type. Use 'none' if no averaging is to be done" << endl;
							throw;
						}	
					} // frames processed...
					
					// aggregate results
					// this corresponds to orientational averaging
					if (ti->mode=="none") {

						// no correlation -> instanenous scattering intensity
						for(map<int,vector<complex<double> > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++) {
							double scatsum=0;				
							for(vector<complex<double> >::iterator si=sbfi->second.begin();si!=sbfi->second.end();si++) {
								scatsum += abs(conj(*si)*(*si));
							}
							qF_Task_results[make_pair(ti->q,sbfi->first)] = scatsum/(sbfi->second.size()); 
						}
					}
					else if (mode=="time") {

						// check size of vector<scattering amplitudes>, this is the size of unfolded q vectors
						size_t aqscount = 0;
						for(map<int,vector<complex<double> > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++) {
							if (aqscount==0) aqscount = sbfi->second.size(); else if (aqscount!=sbfi->second.size()) throw;
						}

						for(size_t asqi = 0; asqi < aqscount; ++asqi)
						{
							// communicate current aqs
							vector<pair<int,complex<double> > > aqs_by_frame;				
							for(map<int,vector<complex<double> > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++)
							{
								aqs_by_frame.push_back(make_pair(sbfi->first,sbfi->second[asqi]));
							}
							vector<vector<pair<int,complex<double> > > > vector_in,vector_out;
							vector_in.resize(local.size(),aqs_by_frame);

							boost::mpi::all_to_all(local,vector_in,vector_out);

							// decompose vector_out into new table: frame <-> Aqs, use frame number as implicit position
							vector<complex<double> > aqs; aqs.resize(frames.size());
							for(size_t i = 0; i < vector_out.size(); ++i)
							{
								for(size_t j = 0; j < vector_out[i].size(); ++j)
								{
									aqs[vector_out[i][j].first]=vector_out[i][j].second;
								}
							}

							// now determine which correlation 'step' WE need to do:
							vector<int> stepsizes = get_step_assignments(local.rank(),local.size(),frames.size()/2); // frames.size()/2 is the length of the correlation we are interested in...

							map<int,complex<double> > my_AAconj;

							for(vector<int>::iterator ssi=stepsizes.begin();ssi!=stepsizes.end();ssi++) {
								complex<double> AAconj_sum=0;
								int AAconj_count=0;					
								size_t current_frame=0;
								while( (current_frame+(*ssi))<frames.size()) {
									AAconj_sum += aqs[current_frame]*conj(aqs[current_frame+(*ssi)]); // do we have to take the "REAL" part only?
									AAconj_count++;
									current_frame++;
								}
								my_AAconj[*ssi] = AAconj_sum / double(AAconj_count);
							}

							vector<map<int,complex<double> > > all_AAconj;
							boost::mpi::gather(local,my_AAconj,all_AAconj,0);

							if (local.rank()==0) {
								// reorder results, use stepsize as implicit index
			//					vector<complex<double> > AAconj_series; AAconj_series.resize((frames.size()/2)+((frames.size()+1)%2)); // correlation length = size

								vector<complex<double> > AAconj_series; AAconj_series.resize((frames.size()/2)+1); // correlation length = size					
								for(vector<map<int,complex<double> > >::iterator ai=all_AAconj.begin();ai!=all_AAconj.end();ai++) {
									for(map<int,complex<double> >::iterator aii=ai->begin();aii!=ai->end();aii++) {
										AAconj_series[aii->first]=aii->second;
									}
								}
			//					ofstream ofile("test.dat",ios_base::app);
								for(size_t i = 0; i < AAconj_series.size(); ++i)
								{
									q_Task_results[make_pair(ti->q,i)] = AAconj_series[i].real(); 
									// check THIS! take the real part or do a conj-multiply?
			//						ofile << ti->q.x << "\t" << ti->q.y << "t" << ti->q.z << "\t" << i << "\t" << AAconj_series[i].real() << endl;
								}
			//					ofile << endl;
							}

						}

						// wait for all nodes to finish, before communicating
						local.barrier();
					}
							
				}
			
//			} // qvectors blocking end

		taskcounter++;
	}
	
	// wait for all nodes to finish their computations....
	world.barrier();

	// make node 0 empty
	if (rank==0) sample.frames.clear_cache();

	// we can use one type of result table right now...
	Scatterdata<double> scat;

	// aggregate all results on node 0
	// distinguish between different modes
	if (mode=="none") {
		vector<map<pair<CartesianCoor3D,size_t>, double > > all_qF_Task_results; // used if mode = "none"
		boost::mpi::gather(world,qF_Task_results,all_qF_Task_results,0);
		// re-sort everything
		for(vector< map<pair<CartesianCoor3D,size_t>, double > >::iterator ar=all_qF_Task_results.begin();ar!=all_qF_Task_results.end();ar++) {
			for(map<pair<CartesianCoor3D,size_t>,double>::iterator mi=ar->begin();mi!=ar->end();mi++) {
				scat[mi->first.first][mi->first.second]=mi->second;				
			}
		}

	} 
	else if (mode=="time") {
		vector< map<pair<CartesianCoor3D,size_t>, double > > all_q_Task_results; // used if mode = "time"
		boost::mpi::gather(world,q_Task_results,all_q_Task_results,0);
		// re-sort everything
		for(vector< map<pair<CartesianCoor3D,size_t>, double > >::iterator ar=all_q_Task_results.begin();ar!=all_q_Task_results.end();ar++) {
			for(map<pair<CartesianCoor3D,size_t>,double>::iterator mi=ar->begin();mi!=ar->end();mi++) {
				scat[mi->first.first][mi->first.second]=mi->second;								
			}
		}
		
	}
	
	//------------------------------------------//
	//
	// Output
	//
	//------------------------------------------//	

	if (rank==0) {
		// make sure after hard work, results are NOT lost...
		string outformat;
		bool outputerror=false;
		try {
		
			if (Settings::get("main").exists("output")) {
				string prefix;
				if (Settings::get("main")["output"].exists("prefix")) {
					prefix = (const char *) Settings::get("main")["output"]["prefix"];
				}
				else {
					prefix = "scattering-data-";
				}
		
				libconfig::Setting* setting = &Settings::get("main")["output"];
				for (int i=0;i<setting->getLength();i++) {
					if ((*setting)[i].getType()!=libconfig::Setting::TypeGroup) continue;
					string settings_name = (*setting)[i].getName();

					if ((*setting)[i].exists("type")) {
						if ((*setting)[i].exists("format")) {
							outformat = (char const *) (*setting)[i]["format"];
						}
						else {
							outformat = "txt";
						}
						stringstream ffname; ffname << prefix << settings_name << "." << outformat;
						string typestring = (char const *) (*setting)[i]["type"];
						clog << "INFO>> " << "Writing results to " << ffname.str() << " via method " << typestring << " in format: " << outformat << endl;

						if (typestring == "plain") {
							scat.write_plain(ffname.str(),outformat);
						}
						if (typestring == "average") {
							scat.write_average(ffname.str(),outformat);
						}
					}
					else {
						outputerror = true;
						continue;
					}
				}		
		
			}
			else {
				outputerror = true;
			}
		
			if (outputerror) {
				clog << "INFO>> " << "Trajectory format not understood. Dumping to standard log." << endl;
				scat.dump(clog.rdbuf());			
			}
		}
		catch (...) {
			cerr << "ERROR>> " << "An Error occured during output. Dumping any results to stdlog" << endl;
			scat.dump(clog.rdbuf());
		}
	
		clog << "INFO>> Successfully finished... Have a nice day!"<< endl;
	
	//------------------------------------------//
	//
	// Send hangups to all compute nodes...
	//
	//------------------------------------------//	
	
//		for (int i=1;i<size;i++) {
//			world.send(i,MPI_TAG_HANGUP);
//		}		
	
	}
	
	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	return 0;
}

// end of file
