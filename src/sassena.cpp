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
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "log.hpp"
#include "parameters.hpp"
#include "sample.hpp"
#include "scatterdata.hpp"
#include "settings.hpp"
#include "tasks.hpp"
#include "timer.hpp"

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
		seconds = (microseconds < 0) ? seconds-1 : seconds;
		microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;		
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

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	Params::Inst();
	
	Info::Inst()->set_prefix(to_s(rank)+string(".Info>>"));
	Warn::Inst()->set_prefix(to_s(rank)+string(".Warn>>"));
	Err::Inst()->set_prefix(to_s(rank)+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();
	
	Timer timer;
	timer.start("total");
	
	if (rank==0) {
		//------------------------------------------//
		//
		// Some welcome message....
		//
		//------------------------------------------//	
	
		Info::Inst()->write("This software is being developed by Benjamin Lindner and Franci Merzel.  ");
		Info::Inst()->write("For help, suggestions or correspondense use:                             ");
		Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
		Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Main Developer, Methodology )            ");
		Info::Inst()->write("For publications use the following references:                           ");
		Info::Inst()->write(".........................................................................");
		Info::Inst()->write("1. SASSIM: a method for calculating small-angle X-ray and                ");
		Info::Inst()->write("   neutron scattering and the associated molecular envelope from         ");
		Info::Inst()->write("   explicit-atom models of solvated proteins, F. Merzel and J. C. Smith, ");
		Info::Inst()->write("   Acta Cryst. (2002). D58, 242-249                                      ");
		Info::Inst()->write(".........................................................................");
		Info::Inst()->write("SASSIM is a software for calculating scattering spectra from all-atomic..");
		Info::Inst()->write("");

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
	
		timer.start("sample::setup");

		params->init(string(argv[1]),"conf");
		
		database->init(string(argv[2]),"conf");

	    // create the sample via structure file	
		Info::Inst()->write(string("Reading structure information from file: ")+params->sample.structure.file);
		sample.add_atoms(params->sample.structure.file,params->sample.structure.format);
		Info::Inst()->write(string("Done. Atoms read: ")+to_s(sample.atoms.size()));

	   	// add selections / groups
		for(map<string,SampleGroupParameters>::iterator sgpi = params->sample.groups.begin();sgpi!=params->sample.groups.end();sgpi++)
		{
			SampleGroupParameters& sp = sgpi->second;
			Info::Inst()->write(string("Reading selection file: ")+sp.file);
			sample.atoms.add_selection(sp.name,sp.file,sp.format,sp.select,sp.select_value);
		}

		// add custom selections here ( if not already set! )
		if (sample.atoms.selections.find("system")==sample.atoms.selections.end()) {
			// this shortcut creates a full selection
			sample.atoms.add_selection("system",true);		
		}
	
		// apply deuteration
		for(size_t i = 0; i < params->sample.deuter.size(); ++i)
		{
			sample.deuter(params->sample.deuter[i]);
		}
			
		// read in frame information
		for(size_t i = 0; i < params->sample.frames.size(); ++i)
		{
			pair<string,string> fnt = params->sample.frames[i];
			Info::Inst()->write(string("Reading frames from: ")+fnt.first);
			size_t nof = sample.frames.add_frameset(fnt.first,fnt.second,sample.atoms);			
			Info::Inst()->write(string("Found ")+to_s(nof)+string(" frames"));			
		}
		Info::Inst()->write(string("Total number of frames found: ")+to_s(sample.frames.size()));
		
		// select wrapping behaviour
		if (params->sample.pbc.wrapping) {
			sample.frames.wrapping=true;
			string centergroup = params->sample.pbc.center;
			sample.atoms.assert_selection(centergroup);
			sample.frames.centergroup_selection = sample.atoms.selections[centergroup];			
		}
		else {
			sample.frames.wrapping = false;
		}

		timer.stop("sample::setup");

		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		timer.start("sample::preparation");

		if (params->scattering.background.method=="auto") {

			// determine background scattering density from grid based analysis
			// necessary:
			// group = name of the group which accounts for the background
			// resolution = grid resolution
			// hydration = hydration layer around each particle, 
			//    ie. solvent molecules within this grid distance are not counted towards bulk
			int r     = params->scattering.background.resolution; // resolution of grid
			double h  = params->scattering.background.hydration; // hydration layer of each particle (not background)
		
			Info::Inst()->write("Calculate background scattering length density");
			Info::Inst()->write(string("using grid resolution ")+to_s(r)+string(" and hydration of ")+to_s(h));
			pair<double,double>	b0 = Analysis::background_avg(sample,r,h,CartesianCoor3D(0,0,0));
			Info::Inst()->write(string("Average background scattering length: ")+to_s(b0.first)+string(" +- ")+to_s(b0.second));
			params->scattering.background.value = b0.first;
			
		}
		sample.background = params->scattering.background.value;
		Info::Inst()->write(string("Set background scattering length density to value of ")+to_s(sample.background));
		
		
		timer.stop("sample::preparation");
		
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

		sample.frames.clear_cache(); // reduce overhead

		Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");

		world.barrier();

		timer.start("sample::communication");

		broadcast(world,*params,0);

		world.barrier();

		broadcast(world,*database,0);

		world.barrier();

		broadcast(world,sample,0);
		
		timer.stop("sample::communication");
		
	}
	else {
		
		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, *params);			

		world.barrier();
	
		world.recv(0,boost::mpi::any_tag, *database);			

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
	std::vector<CartesianCoor3D>& qqqvectors = params->scattering.qqqvectors;
	std::vector<int> frames;
	
	timer.start("task preparation");

	// fill frames
	for(size_t i=0;i<sample.frames.size();i+= params->scattering.framestride) {
		frames.push_back(i);
	}	
	if (rank==0) {
		Info::Inst()->write(string("Framestride used: ")+to_s(params->scattering.framestride));
		Info::Inst()->write(string("Number of frames to be evaluated: ")+to_s(frames.size()));		
	}

	int taskcounter = 0; int progress = 0;
	
	// generate Tasks
	string mode = params->scattering.correlation.type;

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
	string target = params->scattering.target;

	if (rank==0) { 
		Info::Inst()->write(string("Scattering target selection: ")+target);
	}
	
	timer.stop("task preparation");

	gettimeofday(&start, 0);
	int progess = 0;
	
	map<pair<CartesianCoor3D,size_t>, double > qF_Task_results; // used if mode = "none"
	map<pair<CartesianCoor3D,size_t>, double > q_Task_results; // used if mode = "time"
	
	// if mode node
	// this node will aggregrate results for:
	for (Tasks::iterator ti=my_tasks.begin();ti!=my_tasks.end();ti++) {
		vector<int> frames_i = ti->frames(rank);
		for(size_t i = 0; i < frames_i.size(); ++i)
		{
			qF_Task_results[make_pair(ti->q,frames_i[i])]=0;
		}
	}
	// if mode time
	// this node will aggreate results for:
	for (Tasks::iterator ti=my_tasks.begin();ti!=my_tasks.end();ti++) {
		vector<int> stepsizes = get_step_assignments(local.rank(),local.size(),frames.size()/2); // frames.size()/2 is the length of the correlation we are interested in...
		for(size_t i = 0; i < stepsizes.size(); ++i)
		{
			q_Task_results[make_pair(ti->q,stepsizes[i])]=0;
		}
	}
	
//	Info::Inst()->write("result table:") ;
//	for(map<pair<CartesianCoor3D,size_t>, double >::iterator i = q_Task_results.begin(); i != q_Task_results.end(); ++i)
//	{
//		Info::Inst()->write(string("(")+to_s(i->first.first.x)+string(",")+to_s(i->first.first.y)+string(",")+to_s(i->first.first.z)+string("),")+to_s(i->first.second)+string(":")+to_s(i->second)) ;
//	}

	
	
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
			timer.start("scatter");

			timer.start("scatter::scatteramp");

			// set scattering amplitudes for current q vector
			Analysis::set_scatteramp(sample,sample.atoms.selections[target],ti->q,true);

			timer.stop("scatter::scatteramp");
			
			// generate q-vector list for orientational averaging
			vector<CartesianCoor3D> qvectors;
			
			timer.start("scatter::qvector-unfold");
						
			if (params->scattering.average.method=="none") {
				qvectors.push_back(ti->q);
			} 
			else {
			
				double resolution = params->scattering.average.resolution;
				string avtype = params->scattering.average.type;
			
				uint32_t qseed = 1; // ti->qseed			
				if (avtype=="none") {
					qvectors.push_back(ti->q);
				}
				else if (avtype=="sphere") {
					string avmethod = params->scattering.average.method;
					if (avmethod=="bruteforce") {
						
						string avvectors = params->scattering.average.vectors;
						Analysis::qvectors_unfold_sphere(avvectors,ti->q,qseed,resolution,qvectors);
					}
					if (avmethod=="multipole") {
						// don't unfold. 
						// multipole expansion works on one qvector
					}		
				}
				else if (avtype=="cylinder") {
					string avmethod = params->scattering.average.method;
					if (avmethod=="bruteforce") {
						string avvectors = params->scattering.average.vectors;
						Analysis::qvectors_unfold_cylinder(avvectors,ti->q,qseed,resolution,qvectors);
					}
					if (avmethod=="multipole") {
						// don't unfold. 
						// multipole expansion works on one qvector
					}
				}	
			}
			
			timer.stop("scatter::qvector-unfold");

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
				if (params->scattering.interference.type == "self") {

					map<int,vector<vector<complex<double> > > > scatbyframe; // frame -> qvectors/atoms/amplitude

					// make little pieces, otherwise we'll run into memory issues
					vector<Atomselection> target_selections = sample.atoms.selections[target].slice(1);
//					Info::Inst()->write("slicing target selection, sizes:");
//					Info::Inst()->write(string("target_selection.size: ")+to_s(target_selections.size()));
					
//   				stringstream ss;
//   				for(size_t i = 0; i < target_selections.size(); ++i)
//   				{
//   					ss << target_selections[i].size() << " ";
//   				}
//   				Info::Inst()->write(ss.str());

					for(size_t tsi = 0; tsi < target_selections.size(); ++tsi)
					{
						Atomselection& target_selection = target_selections[tsi];
			
						// iterate through all frames this node is supposed to do
						for(size_t i = 0; i < frames_i.size(); ++i)
						{
							
							timer.start("scatter::loadframe");
							
							sample.frames.load(frames_i[i],sample.atoms,target_selection);				
	
							timer.stop("scatter::loadframe");
	
							// holds the scattering amplitudes for the current frame:
							vector<vector<complex<double> > >& scattering_amplitudes = scatbyframe[frames_i[i]];
	
							string avtype = params->scattering.average.type;
	
							timer.start("scatter::compute");
						
							if (avtype=="none") {
								// qseed not used here
								vector<complex<double> > scat;
								Analysis::scatter_none(sample,target_selection,ti->q,scat);
								scattering_amplitudes.push_back(scat);
							}
							else if (avtype=="sphere") {
								string avmethod = params->scattering.average.method;
								if (avmethod=="bruteforce") {
									Analysis::scatter_vectors(sample,target_selection,qvectors,scattering_amplitudes);						
								}
								if (avmethod=="multipole") {
	//								if (resolution==-1.0) resolution=17.0;
									cerr << "ERROR>> " << " Multipole averaging w/ interference: self not (yet) implemented." << endl;
									throw;
	//					   			Analysis::scatter_sphere_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
								}		
							}
							else if (avtype=="cylinder") {
								string avmethod = params->scattering.average.method;	
								if (avmethod=="bruteforce") {
									Analysis::scatter_vectors(sample,target_selection,qvectors,scattering_amplitudes);		
								}
								if (avmethod=="multipole") {
	//								if (resolution==-1.0) resolution=10.0;
									cerr << "ERROR>> " << " Multipole averaging w/ interference: self not (yet) implemented." << endl;							
									throw;
	//					   			Analysis::scatter_cylinder_multipole(sample,as,ti->q,resolution,scattering_amplitudes);			
								}
							}	
							else {
								Err::Inst()->write("unrecognized averaging type. Use 'none' if no averaging is to be done");
								throw;
							}	
	
							timer.stop("scatter::compute");
	
						} // frames processed
	
	
						timer.start("scatter::aggregate");
	
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
								qF_Task_results[make_pair(ti->q,sbfi->first)] += scatsum/(sbfi->second.size()); 
							}
						}
						else if (mode=="time") {
	
							// check size of vector<scattering amplitudes>, this is the size of unfolded q vectors
							size_t aqscount = 0;
							for(map<int,vector<vector<complex<double> > > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++) {
								if (aqscount==0) aqscount = sbfi->second.size(); else if (aqscount!=sbfi->second.size()) throw;
							}
	
							for(size_t asqi = 0; asqi < aqscount; ++asqi) // do 1 unfolded q vector at a time
							{
							// communicate vector of current aqs
	
								// sbfi->first = frame; sbfi->second = qvector/atoms/scattering-amplitudes
								vector<pair<int,vector<complex<double> > > > aq_vectors_by_frame;				
								for(map<int,vector<vector<complex<double> > > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++)
								{
									aq_vectors_by_frame.push_back(make_pair(sbfi->first,sbfi->second[asqi]));
								}
								
								vector<vector<pair<int,vector<complex<double> > > > > vector_out;
	//							vector_in.resize(local.size(),aq_vectors_by_frame);
	
								timer.start("scatter::agg::correlate");
	
								boost::mpi::all_gather(local,aq_vectors_by_frame,vector_out);
	
								timer.start("scatter::agg::corr::comp");
	
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
									
									if (q_Task_results.find(make_pair(ti->q,*ssi))==q_Task_results.end()) throw; // safety switch
//										// check THIS! take the real part or do a conj-multiply?
									q_Task_results[make_pair(ti->q,*ssi)] += AAconj_sum.real() / double(AAconj_count) / double(aqscount);
								}
								
								timer.start("scatter::agg::corr::comp");
								
								timer.start("scatter::agg::correlate");
								
								
							} // this iterates the unfolded q vectors
	
							// wait for all nodes to finish, before communicating
							local.barrier();
						
							timer.stop("scatter::aggregate");
							
						} // this concludes aggregation (mode-if)

					} // this iterates the target selection slices
					
				}
				else if (params->scattering.interference.type == "all") {

					map<int,vector<complex<double> > > scatbyframe; // frame -> qvectors/amplitude
			
					// iterate through all frames this node is supposed to do
					for(size_t i = 0; i < frames_i.size(); ++i)
					{
						timer.start("scatter::loadframe");
											
						sample.frames.load(frames_i[i],sample.atoms,sample.atoms.selections[target]);				

						timer.stop("scatter::loadframe");

//						// introducing the supergrid 3x3x3 elements
//						map< pair<int,pair<int,int> > , complex<double> >  supergrid;
//						int maxn =35;
//						for (int kkk=0;kkk<maxn;kkk++) {
//							for (int lll=0;lll<maxn;lll++) {
//								for (int mmm=0;mmm<maxn;mmm++) {
//									double e1 = kkk * ti->q * sample.frames.current().unitcell[0];
//									double e2 = lll * ti->q * sample.frames.current().unitcell[1];
//									double e3 = mmm * ti->q * sample.frames.current().unitcell[2];
//									int maxi = ( (kkk>lll) ? kkk : lll  ) > mmm ? ( (kkk>lll) ? kkk : lll  ) : mmm;
//									double damping_coeff = -5.0;
//									double damping = exp(damping_coeff * maxi);
//									supergrid[make_pair(kkk,make_pair(lll,mmm))]= exp(complex<double>(0.0,e1+e2+e3)) * damping;
//									
//								}
//							}
//						}
						//

						Atomselection& target_selection = sample.atoms.selections[target];
						// holds the scattering amplitudes for the current frame:
						vector<complex<double> >& scattering_amplitudes = scatbyframe[frames_i[i]];
								
						string avtype = params->scattering.average.type;
						double resolution = params->scattering.average.resolution;
						timer.start("scatter::compute");
												
						if (avtype=="none") {
							// qseed not used here
///							cout << "INFO>> " << "Calculating scattering for: " << ti->q << endl;
							
							complex<double> scat = Analysis::scatter_none(sample,target_selection,ti->q);
							
//							cout << "INFO>> " << "I=: " << scat << endl;
//							complex<double> scatsuperposed(0.0,0.0);
//							for (map< pair<int,pair<int,int> > , complex<double> >::iterator sgi = supergrid.begin(); sgi!=supergrid.end(); sgi++) {
//								scatsuperposed += (scat*sgi->second);
//							}
														
							scattering_amplitudes.push_back(scat);
//							scattering_amplitudes.push_back(scatsuperposed);
						}
						else if (avtype=="sphere") {
							string avmethod = params->scattering.average.method;
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,target_selection,qvectors,scattering_amplitudes);						
							}
							if (avmethod=="multipole") {
					   			Analysis::scatter_sphere_multipole(sample,target_selection,ti->q,resolution,scattering_amplitudes);			
							}		
						}
						else if (avtype=="cylinder") {
							string avmethod = params->scattering.average.method;
							if (avmethod=="bruteforce") {
								Analysis::scatter_vectors(sample,target_selection,qvectors,scattering_amplitudes);		
							}
							if (avmethod=="multipole") {
					   			Analysis::scatter_cylinder_multipole(sample,target_selection,ti->q,resolution,scattering_amplitudes);			
							}
						}	
						else {
							cerr << "ERROR>> " << "unrecognized averaging type. Use 'none' if no averaging is to be done" << endl;
							throw;
						}	
						
						timer.stop("scatter::compute");
						
					} // frames processed...
					
					
					timer.start("scatter::aggregate");
					
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

						//
						int aqsblock = 500;
						for(size_t Ti = 0; Ti < aqscount; Ti+=aqsblock)
						{

						for(size_t asqi = Ti; asqi < ( (aqscount<(Ti+aqsblock-1)) ? aqscount : (Ti+aqsblock-1) ) ; ++asqi)
						{
							// communicate current aqs
							vector<pair<int,complex<double> > > aqs_by_frame;				
							for(map<int,vector<complex<double> > >::iterator sbfi=scatbyframe.begin();sbfi!=scatbyframe.end();sbfi++)
							{
								aqs_by_frame.push_back(make_pair(sbfi->first,sbfi->second[asqi]));
							}
							vector<vector<pair<int,complex<double> > > > vector_out;
//							vector_in.resize(local.size(),aqs_by_frame);

							timer.start("scatter::agg::correlate");

							boost::mpi::all_gather(local,aqs_by_frame,vector_out);

							// decompose vector_out into new table: frame <-> Aqs, use frame number as implicit position
							vector<complex<double> > aqs; aqs.resize(frames.size());
							for(size_t i = 0; i < vector_out.size(); ++i) {
								for(size_t j = 0; j < vector_out[i].size(); ++j)
								{
									aqs[vector_out[i][j].first]=vector_out[i][j].second;
								}
							}

							// now determine which correlation 'step' WE need to do:
							vector<int> stepsizes = get_step_assignments(local.rank(),local.size(),frames.size()/2); // frames.size()/2 is the length of the correlation we are interested in...

							map<int,complex<double> > my_AAconj;

							timer.start("scatter::agg::corr::comp");

							for(vector<int>::iterator ssi=stepsizes.begin();ssi!=stepsizes.end();ssi++) {
								complex<double> AAconj_sum=0;
								int AAconj_count=0;					
								size_t current_frame=0;
								while( (current_frame+(*ssi))<frames.size()) {
									AAconj_sum += aqs[current_frame]*conj(aqs[current_frame+(*ssi)]); // do we have to take the "REAL" part only?
									AAconj_count++;
									current_frame++;
								}
								
								if (q_Task_results.find(make_pair(ti->q,*ssi))==q_Task_results.end()) throw; // safety switch
//										// check THIS! take the real part or do a conj-multiply?
								q_Task_results[make_pair(ti->q,*ssi)] += AAconj_sum.real() / double(AAconj_count) / double(aqscount);
								
							}

							timer.stop("scatter::agg::corr::comp");
							
							timer.stop("scatter::agg::correlate");							

						}
						}
						// wait for all nodes to finish, before communicating
						local.barrier();
						
					}
					timer.stop("scatter::aggregate");
							
				}
			
//			} // qvectors blocking end

		timer.stop("scatter");

		taskcounter++;
	}
	
	// wait for all nodes to finish their computations....
	world.barrier();

	timer.start("result::aggregrate");

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
	
	timer.start("result::aggregrate");
	
	//------------------------------------------//
	//
	// Output
	//
	//------------------------------------------//	

	timer.start("result::output");

	if (rank==0) {
		// make sure after hard work, results are NOT lost...
		string outformat;
		bool outputerror=false;
		try {
			for(size_t i = 0; i < params->output.files.size(); ++i)
			{
				string otype   = params->output.files[i].type;
				string oformat = params->output.files[i].format;
				string ofn = params->output.files[i].filename;					
            
				Info::Inst()->write(string("Writing results to ")+ofn+string(" via method ")+otype+string(" in format: ")+oformat);
				
				if (otype=="plain") {
					scat.write_plain(ofn,oformat);
				}
				if (otype=="average") {
					scat.write_average(ofn,oformat);
				}
				
			}
		}
		catch (...) {
			cerr << "ERROR>> " << "An Error occured during output. Dumping any results to stdlog" << endl;
			scat.dump(clog.rdbuf());
			clog.flush();
		}
	

	//------------------------------------------//
	//
	// Send hangups to all compute nodes...
	//
	//------------------------------------------//	
	
//		for (int i=1;i<size;i++) {
//			world.send(i,MPI_TAG_HANGUP);
//		}		
	
	}
	
	timer.stop("result::output");

	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	timer.stop("total");
	
	if (rank==0) {
		vector<string> keys = timer.keys();
				
		clog << "INFO>> " << "                                                                         " << endl;
		clog << "INFO>> " << "                    Performance Analysis                                 " << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> " << " mean and total runtimes:                                                " << endl;				
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> ";
		clog << setw(31) << " measure |";
		clog << setw(12) << " total |";
		clog << setw(10) << " count |";
		clog << setw(10) << " mean |";
		clog << setw(10) << " stddev ";
		clog << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			clog << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.sum(*ki)            << " |";	
			clog << setiosflags(ios::fixed) << setprecision(0) << setw(8) << timer.count(*ki)          << " |";		
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.mean(*ki)           << " |";
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << sqrt(timer.variance(*ki)) << endl;
		}
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		

		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> " << " watermarks:                                                " << endl;				
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		clog << "INFO>> ";
		clog << setw(31) << " measure |";
		clog << setw(12) << " min |";
		clog << setw(10) << " max ";
		clog << endl;
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			clog << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.min(*ki)            << " |";	
			clog << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.max(*ki)            << endl;
		}
		clog << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
	
	}
	
	if (rank==0) { 	
		clog << "INFO>> Successfully finished... Have a nice day!"<< endl;
	}
	return 0;
}

// end of file
