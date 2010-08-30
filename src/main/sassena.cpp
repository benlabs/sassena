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
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/thread.hpp>


// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include "decomposition/decomposition_plan.hpp"
#include "control.hpp"
#include "log.hpp"
#include "report/performance_analyzer.hpp"
#include "report/progress_reporter.hpp"
#include "report/timer.hpp"
#include "sample/sample.hpp"
#include "scatter_devices/scatter_device_factory.hpp"
#include "scatter_devices/scatter_spectrum.hpp"
#include "threadpool/threadpool.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;


//void monitor_thread(boost::mpi::communicator progress_comm,ScatterDevice* pScatterDevice) {
//    ProgressReporter progress_reporter(progress_comm.size(),0); // this acts like a progress bar
//    double maxprogress = progress_comm.size()*1.0;
//    double totalprogress=0;
//    while (totalprogress!=maxprogress) {
//    	double progress = pScatterDevice->progress();
//    	boost::mpi::all_reduce(progress_comm,progress,totalprogress,std::plus<double>());
//    // progress indicator
//		if (progress_comm.rank()==0) { // only let the world head report...
//			progress_reporter.set(totalprogress);
//			progress_reporter.report(); // only reports if significant progress has been made!
//		}
//        cerr << "monitor" << endl;
//		sleep(2);
//	}
//}

void scatter_thread(boost::mpi::communicator world_comm,ScatterDevice* pScatterDevice) {
	int done = 0;
	int alldone = 0;
    ProgressReporter progress_reporter(world_comm.size(),0); // this acts like a progress bar
	
    double totalprogress =0;
    double progress=0;
    // pre-computing progress. don't enter loop if there is no work to do...
	if (pScatterDevice->progress()==1) done=1;
	boost::mpi::all_reduce(world_comm,done,alldone,std::plus<double>());
	

	if (alldone==world_comm.size()) {
        if (world_comm.rank()==0) {		
		    Info::Inst()->write("No q vectors left to compute");
	    }
        return;
	}
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        if (world_comm.rank()==0) {		
    		Info::Inst()->write("Qvectors orientations used for averaging: ");
            for (int i=0;i<Params::Inst()->scattering.average.orientation.vectors.size();i++) {
                string qvector = "";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].x);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].y);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].z);
    		    Info::Inst()->write(qvector);                                
            }
    	}	    
	}
	
    if (world_comm.rank()==0) {		
		Info::Inst()->write("Starting scattering calculation");
	}
	
	while(alldone!=world_comm.size()) {
		pScatterDevice->compute();
		pScatterDevice->write();
		pScatterDevice->next();
		
    	progress = pScatterDevice->progress();
        boost::mpi::all_reduce(world_comm,progress,totalprogress,std::plus<double>());
        // progress indicator
    	if (world_comm.rank()==0) { // only let the world head report...
    		progress_reporter.set(totalprogress);
    		progress_reporter.report(); // only reports if significant progress has been made!
    	}

		// test total progress
		if (pScatterDevice->progress()==1) done=1;
		boost::mpi::all_reduce(world_comm,done,alldone,std::plus<double>());
	}
}

int main(int argc,char** argv) {

	//------------------------------------------//
	//
	// MPI Initialization
	//
	//------------------------------------------//
  	boost::mpi::environment env(argc, argv);
  	boost::mpi::communicator world;

    // The rank 0 node is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	Params::Inst();
	
	Info::Inst()->set_prefix(to_s(world.rank())+string(".Info>>"));
	Warn::Inst()->set_prefix(to_s(world.rank())+string(".Warn>>"));
	Err::Inst()->set_prefix(to_s(world.rank())+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
	Timer timer;
	timer.start("total");
	
	if (world.rank()==0) {
		//------------------------------------------//
		//
		// Some welcome message....
		//
		//------------------------------------------//	
	
		Info::Inst()->write("This software is being developed by Benjamin Lindner.                    ");
		Info::Inst()->write("For help, suggestions or correspondense use:                             ");
		Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
		Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Methodology)                             ");
		Info::Inst()->write("For publications include the following references:                       ");
		Info::Inst()->write(".........................................................................");
		Info::Inst()->write("1. Sassena - Scattering Calculations on Parallel Computers               ");
		Info::Inst()->write("   to be published                                                       ");		
		Info::Inst()->write(".........................................................................");
        Info::Inst()->write(string("Version Information: ") + string(Sassena_VERSIONSTRING));
		Info::Inst()->write("");
    }

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config", po::value<string>()->default_value("scatter.xml"),  "name of the xml configuration file")
        ("database",po::value<string>()->default_value("db.xml"),  "name of the xml database file")   
        ("output",po::value<string>()->default_value("fqt.h5"),"name of the data output file")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }
    
    if (vm.count("config")==0) {
        Info::Inst()->write("Require configuration file");
        cout << desc << endl;
    }
    if (vm.count("database")==0) {
        Info::Inst()->write("Require database file");
        cout << desc << endl;
    }
    
    bool initstatus = true;
    
	if (world.rank()==0) {
	    
		try {

    		timer.start("sample::setup");
	    
            params->init(vm["config"].as<string>());

            database->init(vm["database"].as<string>());
	    
            sample.init();
		
		    timer.stop("sample::setup");
		    
        } catch (boost::exception const& e ) {
            initstatus = false; 
            Err::Inst()->write("Caught BOOST error, sending hangup to all nodes");
            stringstream ss; ss << diagnostic_information(e);
            Err::Inst()->write(string("Diagnotic information: ") + ss.str());
        } catch (std::exception const& e) {
            initstatus = false; 
            Err::Inst()->write("Caught STD error, sending hangup to all nodes");
            Err::Inst()->write(string("what() : ") + e.what());
        } catch (...) {
            initstatus = false; 
            Err::Inst()->write("Caught error: UNKNOWN sending hangup to all nodes");
        }
        
		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
    }
    
    broadcast(world,initstatus,0);            	
	
	// if something went wrong during initialization, exit now.
    if (!initstatus) return 1;
	
	if (world.rank()==0) {
	    
		Info::Inst()->write(string("Set background scattering length density set to ")+to_s(Params::Inst()->scattering.background.factor));
		
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

	    Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
    }
    
	world.barrier();

	timer.start("sample::communication");

	broadcast(world,*params,0);

	world.barrier();

	broadcast(world,*database,0);

	world.barrier();

	broadcast(world,sample,0);

	world.barrier();
	
	timer.stop("sample::communication");


	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//

//    boost::mpi::communicator progress_comm = boost::mpi::communicator(world,world.group());
    boost::mpi::communicator mutex_comm = boost::mpi::communicator(world,world.group());

    boost::mpi::communicator scatter_comm = boost::mpi::communicator(world,world.group());

	world.barrier();

    if (world.rank()==0)
	    Info::Inst()->write("Setting up parallel environment...");		    
    ScatterDevice* p_ScatterDevice =NULL;
    p_ScatterDevice = ScatterDeviceFactory::create(
    		scatter_comm,
    		sample,
    		params->scattering.qvectors,
    		vm["output"].as<string>());

	world.barrier();

    // start computation tasks
    if (world.rank()==0)
		Info::Inst()->write("Starting computation...");		    
//    boost::thread pthread(boost::bind(&monitor_thread,progress_comm,p_ScatterDevice));
    boost::thread sthread(boost::bind(&scatter_thread,mutex_comm,p_ScatterDevice));

    // do a wait on the threads here
//    pthread.join();
    sthread.join();

	if (world.rank()==0) {
		Info::Inst()->write(string("Waiting for all nodes to finish calculations..."));		
	}
    world.barrier();
	
	
	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}
	
	timer.stop("total");
	
    PerformanceAnalyzer perfanal(world,p_ScatterDevice->timer); // collect timing information from everybody.
    delete p_ScatterDevice;
    if (world.rank()==0) {
		perfanal.report();
		perfanal.report_relative(timer.sum("total")*world.size());
    	Info::Inst()->write(string("Total runtime (s): ")+to_s(timer.sum("total")));
    	Info::Inst()->write("Successfully finished... Have a nice day!");
	}




	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	

	return 0;
}


// end of file
