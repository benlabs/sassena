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
      
	// MPI communication is implemented via a master-slave technique.
	// The head node is responsible for the delegation of computations
    // and the collection of data. 
    // The head node also is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:
    // cerr << "ERROR>> (" << rank << "/" << size ")" << " some error message " << endl;
	// clog << "INFO>> " << "############ This is a debug hello from node ############             " << rank << endl;
	if (rank!=0) {
		// clog << "INFO>> " << "############ Node " << rank << " is compute node       "  << endl;		
		main_slave(env,world,argc,argv);
		return 0;
	}


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

	if (size<2) {
		clog << "ERROR>> " << "This version is mpi-enabled. I need at least 2 nodes (master/slave) " << endl;
		clog << "ERROR>> " << "Try again with mpirun -n 2 sassim ...                               " << endl;		
		return -1;
	} 
	
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
   Sample sample( Settings::get_filepath(Settings::get("main")["sample"]["structure"]["file"]), Settings::get("main")["sample"]["structure"]["format"]);

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
		sample.add_selection("system",sample.atoms.begin(),sample.atoms.end());		
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
		sample.add_frame(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"]);
	}

	// select wrapping behaviour
	if (Settings::get("main")["sample"]["pbc"]["wrapping"]) {
		sample.wrapping = true;
		sample.centergroup = (const char *) Settings::get("main")["sample"]["pbc"]["center"];
		clog << "INFO>> " << "Turned wrapping ON with center group " << sample.centergroup << endl;				
	}
	else {
		sample.wrapping = false;		
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

	clog << "INFO>> " << "Exchanging sample information with compute nodes... " << endl;
	broadcast(world,sample,0);
	

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
	
	Scatterdata<double> scat;
	
	// fill qqqvectors
	Settings::get_qqqvectors(qqqvectors);
	
	// fill frames
	int framestride = 1;
	if (Settings::get("main")["scattering"].exists("framestride")) {
		framestride = Settings::get("main")["scattering"]["framestride"];
	}
	for(int i=0;i<sample.dcdframes.number_of_frames();i+= framestride) {
		frames.push_back(i);
	}	

	int currentjobs = 0;
	
	int numberofjobs = qqqvectors.size()*frames.size();
	int jobcounter = 0; int progress = 0;
	clog << "INFO>> " << "Total number of jobs: " << numberofjobs << " (qqqvectors: " << qqqvectors.size() << ", frames: " << frames.size() << ")"<< endl;
//	clog << "INFO>> " << "Progress(%): ";
	gettimeofday(&start, 0);
	// this logic provides approx. 20 submissions for each node  
	// however, the maximum number of simultaneous keys is limited by the second layer of the loop: qqqvector currently. 
	// this is a performance variable.... subject to improvement
	// it's a trade off: load balancing vs. stress for the master loop
	// indicator for bad master loop performance: currentjobs << size-1
	// extreme case: numsimjobs = numberofjobs/(size-1)  , unless not limited by 2nd loop provides NO LOAD BALANCING
	// factor 0.05 = 1/20 provides load balancing with a granularity of 20. i.e. it can handle a bad case scenario of 2000% load imbalance (overestimated)
	// -> for 160 cpus, 1000 qqqvectors, 50 frames. 0.05*50000/159 = 15.7
	// -> for 160 cpus, 1000 qqqvectors, 5000 frames. 0.05*5000000/159 = 1570 -- 1000 limit!	
	// -> for 160 cpus, 1000000 qqqvectors, 1 frames. 0.05*1000000/159 = 314	
	int numsimjobs = int(0.05*numberofjobs/(size-1));	
	numsimjobs = (numsimjobs>1) ? numsimjobs : 1 ; 
	
	for(vector<int>::iterator fi=frames.begin();fi!=frames.end();fi++) {
	
		for (vector<CartesianCoor3D>::iterator qi=qqqvectors.begin();qi!=qqqvectors.end();) {
		

			int percent = 	int(jobcounter*100/numberofjobs)  ;	
			if ( percent >= progress  ) {
				gettimeofday(&end, 0);				
				int seconds = (end.tv_sec-start.tv_sec);
				int microseconds = (end.tv_usec-start.tv_usec);
				microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
				seconds = (microseconds < 0) ? seconds-1 : seconds;								
				double etotal = 100 * ( seconds + (microseconds/1000000.0) ) / percent;
				double eta = etotal - seconds + (microseconds/1000000.0);
				clog << "INFO>> " << "Progress: " << percent << "%" << " , ETOTAL(s): " <<  etotal  << " , ETA(s): " <<  eta  << endl;
				clog << "INFO>> " << "Current load (currentjobs): " <<  currentjobs  << endl;			
				if (percent<10) {
					progress = 1*(percent/1) + 1;									
				}
				else {
					progress = 10*(percent/10) + 10;									
				}				
			}
			
			boost::mpi::status s;
			bool jobsubmitted = false;	
			while ( !jobsubmitted ) {
				s = world.recv(boost::mpi::any_source, boost::mpi::any_tag);
				
				if (s.tag()==MPI_TAG_IDLE) {
						vector<ScatterdataKey> keys;
						for (int i=0;i<numsimjobs;i++) {
							if (qi==qqqvectors.end()) break;
							keys.push_back(ScatterdataKey(*qi,*fi)); 
							qi++;
						}	
						world.send(s.source(),MPI_TAG_SCATTER);				
						world.recv(s.source(),MPI_TAG_READY);													
						world.send(s.source(),MPI_TAG_SCATTER,keys);							
//						clog << "INFO>> " << "Sending " << keys.size() << " jobs to node " << s.source() << endl;											
//						if (keys.size()!=numsimjobs) clog << "INFO>> " << "Maximum jobs per node limited" << s.source() << endl;																	
						// increase jobcounter
						// send job to compute node
						currentjobs++;
						jobcounter+=keys.size();
						jobsubmitted=true;
				}
				
				else if (s.tag()==MPI_TAG_SCATTERRESULT) {
					world.send(s.source(),MPI_TAG_READY);
					vector<pair<ScatterdataKey,double> > keyvals;				
					world.recv(s.source(),MPI_TAG_SCATTERRESULT,keyvals);
//					clog << "INFO>> " << "Received " << keyvals.size() << " results from node " << s.source() << endl;					
					for (vector<pair<ScatterdataKey,double> >::iterator kvi=keyvals.begin();kvi!=keyvals.end();kvi++) {
						scat[kvi->first]=kvi->second;
					}				
					currentjobs--;
				} 
				else {
					// dunno what to do....
					// get next message
				}
			}


						// this one has to be done by the compute node!
			//			pair<CartesianCoor3D,int> key(*qi,i);	
			//			sample.read_frame(i);
			//			Analysis::set_scatteramp(sample,sample.atomselections[target],*qi,true,probe);
			//			pair<CartesianCoor3D,int> key(*qi,i);
			//			scat[key]=Analysis::scatter(sample,sample.atomselections[target],*qi);


		}
		/////////////////////			
		// END Frame loop part
		/////////////////////		
	}

		
	//catch pending results...	
	while (currentjobs>0) {
		boost::mpi::status s;
		s = world.recv(boost::mpi::any_source, boost::mpi::any_tag);
		
		if (s.tag()==MPI_TAG_SCATTERRESULT) {
			world.send(s.source(),MPI_TAG_READY);
			vector<pair<ScatterdataKey,double> > keyvals;				
			world.recv(s.source(),MPI_TAG_SCATTERRESULT,keyvals);
//			clog << "INFO>> " << "Received " << keyvals.size() << " results from node " << s.source() << endl;
			for (vector<pair<ScatterdataKey,double> >::iterator kvi=keyvals.begin();kvi!=keyvals.end();kvi++) {
				scat[kvi->first]=kvi->second;
			}				
			currentjobs--;
		}	
	}
	
	//------------------------------------------//
	//
	// Output
	//
	//------------------------------------------//	

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
	
	for (int i=1;i<size;i++) {
		world.send(i,MPI_TAG_HANGUP);
	}		
	
	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	return 0;
}

void main_slave(boost::mpi::environment& env, boost::mpi::communicator& world,int argc,char** argv) {

	int rank = world.rank(); int size = world.size();

	// for timing:
 	timeval start, end;	


	// read in configuration file
	if (!Settings::read(argc,argv)) { 
		cerr << "ERROR>> (" << rank << "/" << size << ")" << " Error reading the configuration" << endl;
		return;
	}
	
	// create an empty sample
	Sample sample;
	
//	string msg;
//	clog << "INFO>> (" << rank << "/" << size << ")" << "Wait to receive msg: " << msg << endl;	
	world.recv(0,boost::mpi::any_tag, sample);
	
//	ofstream ofs("test-compute.txt");
//	boost::archive::text_oarchive oa(ofs);
//	oa & sample;
	
	
	// implement this as an infinite loop
	// wait for a signal from the head to jump out of it
	while (true) {		
		world.send(0, MPI_TAG_IDLE);
		boost::mpi::status s = world.recv(0,boost::mpi::any_tag);	
			
		if (s.tag()==MPI_TAG_HANGUP) {	
			break;
		}		
		else if (s.tag()==MPI_TAG_SCATTER) {
			boost::mpi::status s2;
			world.send(0, MPI_TAG_READY);
			
			vector<ScatterdataKey> keys;
			vector<pair<ScatterdataKey,double> > keyvals;		
			world.recv(0,MPI_TAG_SCATTER,keys);
			string target = Settings::get("main")["scattering"]["target"];
			
			for (vector<ScatterdataKey>::iterator ki=keys.begin();ki!=keys.end();ki++) {
							
				sample.read_frame(ki->second);				
				gettimeofday(&start, 0);
				Analysis::set_scatteramp(sample,sample.atomselections[target],ki->first,true);
				
				/* debug code for 'rotating the scattering object in situ */
				
//				for (Atomselection::iterator asi=sample.atomselections[target].begin();asi!=sample.atomselections[target].end();asi++) {
//					SphericalCoor3D c2 = sample.currentframe().coord3D(*asi);
//					CartesianCoor3D c = SphericalCoor3D(c2.r,c2.phi,c2.theta+0.25*M_PI);
//					CartesianCoor3D r = c2;
//					r = rotate(r,"y",0.25*M_PI);
//					r = rotate(r,"z",0.25*M_PI);
//					r = rotate(r,"x",0.25*M_PI);					
					
//					CartesianCoor3D c = r;				
//					CartesianCoor3D c = rotate(rotate(rotate(c2,"z",0.1*M_PI),"y",0.1*M_PI),"x",0.1*M_PI);
//					sample.currentframe().x[*asi] = CartesianCoor3D(c).x;
//					sample.currentframe().y[*asi] = CartesianCoor3D(c).y;
//					sample.currentframe().z[*asi] = CartesianCoor3D(c).z;
//				}
//				if (rank==100) {
//					stringstream fn; 
//					fn << "test-" << rank << ".pdb";
//					sample.atoms.write(fn.str(),sample.currentframe());
//				}
				
				/* end of debug code */ 
				
				
				gettimeofday(&end, 0);
				int seconds = (end.tv_sec-start.tv_sec);
				int microseconds = (end.tv_usec-start.tv_usec);
				microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
				seconds = (microseconds < 0) ? seconds-1 : seconds;								
				double tdif = seconds + (microseconds/1000000.0);
	//			clog << "INFO>> " << "Timing set_scatteramp. Time difference: " << tdif  << endl;
				gettimeofday(&start, 0);			
				double scat =Analysis::scatter(sample,sample.atomselections[target],ki->first);		
				gettimeofday(&end, 0);
			
				 seconds = (end.tv_sec-start.tv_sec);
				 microseconds = (end.tv_usec-start.tv_usec);
				microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
				seconds = (microseconds < 0) ? seconds-1 : seconds;								
				 tdif = seconds + (microseconds/1000000.0);
	//			clog << "INFO>> " << "Timing scatter. Time difference: " << tdif  << endl;

				keyvals.push_back(make_pair(*ki,scat));
			}			
			// finished:
			world.send(0, MPI_TAG_SCATTERRESULT);
			world.recv(0, MPI_TAG_READY);			
			world.send(0, MPI_TAG_SCATTERRESULT, keyvals);
		} 
		
		// MPI_RECEIVE
		// if abort -> break
		// if compute -> compute w/ parameters, then send back: results ready, MPI_RECEIVE and send results, repeat loop
		
	}

}

// end of file
