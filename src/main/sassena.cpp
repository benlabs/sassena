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
#include "scatter_devices/io/h5_fqt_interface.hpp"
#include "threadpool/threadpool.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;


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


    vector<size_t> qindexes;
	vector<size_t> frames;
		
	for(size_t i = 0; i < sample.coordinate_sets.size(); ++i)
	{
		frames.push_back(i);
	}
    
    if (world.rank()==0) {
        qindexes = H5FQTInterface::init(vm["output"].as<string>(),params->scattering.qvectors,frames.size());
        
        size_t nq = qindexes.size();
        broadcast(world,nq,0);

        broadcast(world,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);
        
        world.barrier();        
    } else {
        size_t nq=0;
        broadcast(world,nq,0);
        qindexes.resize(nq);
        broadcast(world,reinterpret_cast<size_t*>(&qindexes[0]),qindexes.size(),0);
        world.barrier();        
    }
        
	if (world.rank()==0) {	
    	Info::Inst()->write(string("Searching for decomposition plan: "));
	    Info::Inst()->write(string("nodes    = ")+ to_s(world.size()));
	    Info::Inst()->write(string("qvectors = ")+ to_s(qindexes.size()));
	    Info::Inst()->write(string("frames   = ")+ to_s(sample.coordinate_sets.size()));	
    }
    
	DecompositionPlan dplan(world,qindexes,frames);
	if (world.rank()==0) {	
		Info::Inst()->write(string("Decomposition has ")+to_s(dplan.worlds())+string(" partitions"));
		Info::Inst()->write(string("Static imbalance factor (0 is best): ")+ to_s(dplan.static_imbalance()));
	}
	
	boost::mpi::communicator local = dplan.split();
	
	vector<size_t> myqindexes = dplan.qindexes();

	ScatterDevice* p_ScatterDevice = ScatterDeviceFactory::create(local,sample);

	ScatterSpectrum scatter_spectrum;

    vector<std::complex<double> > fqt;

	ProgressReporter progress_reporter(myqindexes.size(),0); // this acts like a progress bar

    std::vector<CartesianCoor3D> localqvectors;
    
    if (local.rank()==0) 
        localqvectors = H5FQTInterface::get_qvectors(vm["output"].as<string>(),myqindexes);


    world.barrier();
    
    // create communicator for local heads
    // enables file locking
    size_t color = 0;
    if (local.rank()==0) color=1;
    if (myqindexes.size()==0) color=0; // disregard the ones which don't participate in the computation
    boost::mpi::communicator mutexcomm = world.split(color);
	
    std::queue<std::pair<size_t,std::vector<std::complex<double> > > > fqt_queue;
    
	for(size_t qi = 0; qi < myqindexes.size(); ++qi)
	{			
			// progress indicator
			if (world.rank()==0) { // only let the world head report...
				progress_reporter.set(qi);
				progress_reporter.report(); // only reports if significant progress has been made!
			}
			

            CartesianCoor3D qvector;
            if (local.rank()==0) qvector = localqvectors[qi];
            broadcast(local,reinterpret_cast<double*>(&qvector),3,0);

			timer.start("scatter");
			p_ScatterDevice->execute(qvector);
			timer.stop("scatter");
            local.barrier();
			
			if (local.rank()==0) {
    			for(int i = 0; i < mutexcomm.size(); ++i)
    			{
    		        if (mutexcomm.rank()==i) {
    		             H5FQTInterface::store(vm["output"].as<string>(),myqindexes[qi],p_ScatterDevice->get_spectrum());
    		        }
                    mutexcomm.barrier();
    			}
			}
            
            local.barrier();
	}

	if (world.rank()==0) {
		Info::Inst()->write(string("Waiting for all nodes to finish calculations..."));		
	}

    world.barrier();
	
	
	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}

	PerformanceAnalyzer perfanal(world,p_ScatterDevice->timer); // collect timing information from everybody.
	
	delete p_ScatterDevice;

	//------------------------------------------//
	//
	// Finished
	//
	//------------------------------------------//	
	
	timer.stop("total");
		
	if (world.rank()==0) {
		perfanal.report();
		perfanal.report_relative(timer.sum("total")*world.size());
		
		Info::Inst()->write(string("Total runtime (s): ")+to_s(timer.sum("total")));
		Info::Inst()->write("Successfully finished... Have a nice day!");
	}

	return 0;
}

// end of file
