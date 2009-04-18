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

int main(int argc,char** argv) {

	//------------------------------------------//
	//
	// Setting up the sample
	//
	//------------------------------------------//	
	
	clog << "INFO>> " << "Reading configuration file " << argv[1] << endl;
	// read in configuration file, delete command line arguments
	if (!Settings::read(argc,argv)) { cerr << "Error reading the configuration" << endl; throw; }

#ifdef NDEBUG
	clog << "INFO>> " << "NON-DEBUG MODE ACTIVATED" << endl;
#else 
	clog << "INFO>> " << "DEBUG MODE ACTIVATED" << endl;
#endif

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
	
	string target = Settings::get("main")["scattering"]["target"];
	Atomselection& target_sel = sample.atomselections[target];
	clog << "INFO>> " << "compute scattering for " << target_sel.size() << " particles" <<endl;
	
	gettimeofday(&start, 0);
				
	for(vector<int>::iterator fi=frames.begin();fi!=frames.end();fi++) {
		sample.read_frame(*fi);	
		for (vector<CartesianCoor3D>::iterator qi=qqqvectors.begin();qi!=qqqvectors.end();) {
			clog << "INFO>> " << "computing q=" << *qi <<endl;
			Analysis::set_scatteramp(sample,target_sel,*qi,true);
			double scatv =Analysis::scatter(sample,target_sel,*qi);		
			scat[ScatterdataKey(*qi,*fi)]=scatv;
			qi++;
		}
		/////////////////////			
		// END Frame loop part
		/////////////////////		
	}
	
	gettimeofday(&end, 0);				
	int seconds = (end.tv_sec-start.tv_sec);
	int microseconds = (end.tv_usec-start.tv_usec);
	microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
	seconds = (microseconds < 0) ? seconds-1 : seconds;								
	clog << "INFO>> " << "Total time for scattering: " << seconds << "." << microseconds << endl;

	
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
	// Finished
	//
	//------------------------------------------//	
	
	return 0;
}