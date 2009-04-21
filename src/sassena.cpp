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
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <libconfig.h++>
#include <boost/math/special_functions.hpp>

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
	        sample.frames.add_frameset(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"],sample.atoms);
	//		ifstream dcdfilelist(Settings::get_filepath(dcdfilename).c_str());
	//		std::string line;
	//
	//		while (getline(dcdfilelist,line)) {
	//			add_file(line,atoms,(recursion_trigger-1));
	//		}
	//		sample.add_frame(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"]);
		}


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
	string mode = "none";
	Tasks tasks(frames,qqqvectors,world.size(),mode);
	if (rank==0) { 
		clog << "INFO>> " << "Task assignment:" << endl;		
		tasks.print(); 
	}
	
	Tasks my_tasks = tasks.scope(rank);
		
	vector<pair<ScatterdataKey,double> > keyvals;
	string target = Settings::get("main")["scattering"]["target"];

	gettimeofday(&start, 0);
	
	for (Tasks::iterator ti=my_tasks.begin();ti!=my_tasks.end();ti++) {

		if (rank==0) {
			int percent = int(taskcounter*100/my_tasks.size())  ;	
			if ( percent >= progress  ) {
				gettimeofday(&end, 0);				
				int seconds = (end.tv_sec-start.tv_sec);
				int microseconds = (end.tv_usec-start.tv_usec);
				microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
				seconds = (microseconds < 0) ? seconds-1 : seconds;								
				double etotal = 100 * ( seconds + (microseconds/1000000.0) ) / percent;
				double eta = etotal - seconds + (microseconds/1000000.0);
				clog << "INFO>> " << "Progress: " << percent << "%" << " , ETOTAL(s): " <<  etotal  << " , ETA(s): " <<  eta  << endl;			
				if (percent<10) {
					progress = 1*(percent/1) + 1;									
				}
				else {
					progress = 10*(percent/10) + 10;									
				}				
			}	
		}

		
		int frame = ti->rftable[0].second;
		sample.frames.load(frame,sample.atoms,sample.atomselections[target]);				
		Analysis::set_scatteramp(sample,sample.atomselections[target],ti->q,true);

		if (ti->mode=="none") {
			vector<complex<double> > scat;
			// in "none" mode, the scattering intensity is simply the squared amplitude
			uint32_t qseed = 1; // ti->qseed
			Analysis::scatter(sample,sample.atomselections[target],ti->q,qseed,scat);		
			
			double scatsum=0;
			for(size_t i = 0; i < scat.size(); ++i)
			{
				scatsum += abs(conj(scat[i])*scat[i]);
			}
			keyvals.push_back(make_pair(ScatterdataKey(ti->q,frame),scatsum));
			
		}
		else {
			// holds the scattering amplitudes for the current frame:
			vector<complex<double> > scat;
			// when doing dynamic scattering
			// each task has associated frames
			// these define the neighborhood
			// a commuincation scheme has to aggregrate per-vector results on one node
			// then do a FFT convolution
			//the scattering amplitudes have to be aggregrated, then communicated
//			Analysis::scatter(sample,sample.atomselections[target],ti->q,t->qseed,scat);		
//			
//			double scatsum=0;
//			for(size_t i = 0; i < scat.size(); ++i)
//			{
//				scatsum += abs(conj(scat[i])*scat[i]);
//			}
//			keyvals.push_back(make_pair(ScatterdataKey(qqq,frame),scatsum));
			
			
		}

		taskcounter++;
	}
	
	// wait for all nodes to finish their computations....
	world.barrier();
	
	Scatterdata<double> scat;
	
	// aggregate results on node 0
	if (rank==0) {
		for (vector<pair<ScatterdataKey,double> >::iterator kvi=keyvals.begin();kvi!=keyvals.end();kvi++) {
			scat[kvi->first]=kvi->second;
		}		
		for(int i = 1; i < size; ++i)
		{
			// receive keyvals and push to own array
			
			vector<pair<ScatterdataKey,double> > local_keyvals;				
			world.recv(i,MPI_TAG_SCATTERRESULT,local_keyvals);
			for (vector<pair<ScatterdataKey,double> >::iterator kvi=local_keyvals.begin();kvi!=local_keyvals.end();kvi++) {
				scat[kvi->first]=kvi->second;
			}				
			
		}
	}
	else {
		world.send(0, MPI_TAG_SCATTERRESULT, keyvals);
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
				clog << "INFO>> " << "FORMAT: qx qy qz frame I" << endl;		
				scat.dump(clog.rdbuf());			
			}
		}
		catch (...) {
			cerr << "ERROR>> " << "An Error occured during output. Dumping any results to stdlog" << endl;
			clog << "INFO>> " << "FORMAT: qx qy qz frame I" << endl;		
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
