/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
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
#include <boost/asio.hpp>
#include <boost/date_time.hpp>
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
#include "report/timer.hpp"
#include "sample/sample.hpp"
#include "scatter_devices/scatter_device_factory.hpp"
#include "services.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;

void print_title() {
    
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
void print_description() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("......................D.E.S.C.R.I.P.T.I.O.N......................");
	Info::Inst()->write(".................................................................");	

	Info::Inst()->write("This binary computes the scattering intensities directly from");	
	Info::Inst()->write("a molecular dynamics trajectory. ");	
	Info::Inst()->write(".................................................................");	
}



void print_initialization() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("...................I.N.I.T.I.A.L.I.Z.A.T.I.O.N...................");
	Info::Inst()->write(".................................................................");	
}

void print_computation() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write("......................C.O.M.P.U.T.A.T.I.O.N......................");
	Info::Inst()->write(".................................................................");	
}

void print_analysis() {
    Info::Inst()->write(".................................................................");
	Info::Inst()->write(".....................R.U.N...A.N.A.L.Y.S.I.S.....................");
	Info::Inst()->write(".................................................................");	
}



bool init_commandline(int argc,char** argv,po::variables_map& vm) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config", po::value<string>()->default_value("scatter.xml"),  "name of the xml configuration file")
        ("database",po::value<string>()->default_value("db.xml"),  "name of the xml database file")   

        ("scattering-signal-file",po::value<string>()->default_value("signal.h5"),"name of the data output file")
//        ("limits-computation-threads",po::value<string>()->default_value("1"),"threadpool size for orientational averaging")        
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    
    if (vm.find("help")!=vm.end()) {
        cout << desc << endl;
        return false;
    }

    if (vm["config"].defaulted()) {
        Info::Inst()->write("No configuration file specified. Will try to read from scatter.xml");
    }
    if (!boost::filesystem::exists(vm["config"].as<string>())) {
        Err::Inst()->write(vm["config"].as<string>()+string(" does not exist!"));            
        return false;
    }
        
    if (vm["database"].defaulted()) {
        Info::Inst()->write("No database file specified. Will try to read from db.xml");
    }
    if (!boost::filesystem::exists(vm["database"].as<string>())) {
        Err::Inst()->write(vm["database"].as<string>()+string(" does not exist!"));            
        return false;
    }
    
    return true;
}

void read_parameters(po::variables_map vm) {
    Info::Inst()->write("Checking command line for parameter overwrite");
    // second stage , allow command line parameters to overwrite defaults in the configuration file
//	if (!vm["limits-computation-threads"].defaulted()) {
//        Info::Inst()->write(string("Overwrite: limits.computation.threads=")+vm["limits-computation-threads"].as<string>());
 //       Params::Inst()->limits.computation.threads = boost::lexical_cast<size_t>(vm["limits-computation-threads"].as<string>());
   // }
    if (!vm["scattering-signal-file"].defaulted()) {
        Info::Inst()->write(string("Overwrite: scattering.signal.file=")+vm["scattering-signal-file"].as<string>());
        Params::Inst()->scattering.signal.file = vm["scattering-signal-file"].as<string>();        
    }
}

void detect_parameters() {
    if (boost::thread::hardware_concurrency()!=0) {
        Info::Inst()->write(string("Detect: Number of Processors per machine: ")+boost::lexical_cast<string>(boost::thread::hardware_concurrency()));   
    }
}

int main(int argc,char* argv[]) {

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
	
	Info::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Info>>"));
	Warn::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Warn>>"));
	Err::Inst()->set_prefix(boost::lexical_cast<string>(world.rank())+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
	Timer timer;
	timer.start("total");

    bool initstatus = true;
		
    po::variables_map vm;    
    if (world.rank()==0) {
        print_title();
        print_description();
        initstatus = init_commandline(argc,argv,vm);
    }


    broadcast(world,initstatus,0);            	
	// if something went wrong during initialization, exit now.
    if (!initstatus) return 1;

    if (world.rank()==0) print_initialization();
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
    }
    
    broadcast(world,initstatus,0);            	
	// if something went wrong during initialization, exit now.
    if (!initstatus) return 1;
    
    if (world.rank()==0) detect_parameters();
    
    if (world.rank()==0) read_parameters(vm);
    
	if (world.rank()==0) Info::Inst()->write(string("Set background scattering length density set to ")+boost::lexical_cast<string>(Params::Inst()->scattering.background.factor));
		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

	if (world.rank()==0) Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
    
	world.barrier();

	timer.start("sample::communication");
    if (world.rank()==0) Info::Inst()->write("params... ");

	broadcast(world,*params,0);

	world.barrier();

    if (world.rank()==0) Info::Inst()->write("database... ");

	broadcast(world,*database,0);

	world.barrier();
    if (world.rank()==0) Info::Inst()->write("sample... ");

	broadcast(world,sample,0);

	world.barrier();
	
	timer.stop("sample::communication");


	//------------------------------------------//
	//
	// Scattering calculation
	//
	//------------------------------------------//

    if (world.rank()==0) print_computation();

    boost::mpi::communicator scatter_comm = boost::mpi::communicator(world,world.group());

	world.barrier();

    if (world.rank()==0) Info::Inst()->write("Setting up services...");	
    
    boost::asio::io_service io_service;
    
    // setup hdf5writer service and broadcast service to clients
    HDF5WriterService* p_hdf5writer = NULL;    
    MonitorService* p_monitorservice = NULL;    

    // prepare services
    if (world.rank()==0) {
        Info::Inst()->write("Initializing data file service...");	
        p_hdf5writer = new HDF5WriterService(io_service, Params::Inst()->scattering.signal.file,sample.coordinate_sets.size());
        Info::Inst()->write("Initializing monitor service...");	
        p_monitorservice = new MonitorService(io_service, 0, scatter_comm.size());        
    }
    
    // start services
    if (world.rank()==0) {
        Info::Inst()->write("Starting data file service...");	
        p_hdf5writer->run();
        Info::Inst()->write("Starting monitor service...");	
        p_monitorservice->run();
        Info::Inst()->write("Services setup and running...");	
    }

    if (world.rank()==0) {
    }
    // only compute those q vectors which have not been written so far
    // shutdown gracefully if no work is to be done
    initstatus = true;
    vector<CartesianCoor3D> qvectors;
    if (world.rank()==0) {
        
        Info::Inst()->write(string("Checking ")+Params::Inst()->scattering.signal.file + string(" for old results."));
        vector<CartesianCoor3D> oldqvectors = p_hdf5writer->get_qvectors();
        set<CartesianCoor3D> oldqvectors_set;
        for(size_t i = 0; i < oldqvectors.size(); ++i)
        {
            oldqvectors_set.insert(oldqvectors[i]);
        }

        for(size_t i = 0; i < Params::Inst()->scattering.qvectors.size(); ++i)
        {
            if (oldqvectors_set.find(Params::Inst()->scattering.qvectors[i])==oldqvectors_set.end()) {
                qvectors.push_back(Params::Inst()->scattering.qvectors[i]);
            }
        }
        
        if (qvectors.size()!=Params::Inst()->scattering.qvectors.size()) {
            Info::Inst()->write(string("Found ")+ boost::lexical_cast<string>(Params::Inst()->scattering.qvectors.size()-qvectors.size()) + string(" old qvectors."));            
            Info::Inst()->write(string("Total number of qvectors to compute reduced from ")+boost::lexical_cast<string>(Params::Inst()->scattering.qvectors.size()) + string(" to ") + boost::lexical_cast<string>(qvectors.size()) + string("."));            
        }

        if (qvectors.size()==0) {
            Info::Inst()->write("Very good. No qvectors left to compute.");

            Info::Inst()->write("Stopping services.");
            p_monitorservice->hangup();
            p_hdf5writer->hangup();

            initstatus = false;
            Info::Inst()->write("Broadcasting hangup to all nodes. Have a nice day!");
        }
    }
    broadcast(world,initstatus,0);            	
    if (!initstatus) return 1;
    
    // communicate some global parameters
    size_t nq=qvectors.size();
    broadcast(scatter_comm,nq,0);
    if (scatter_comm.rank()!=0) {
        qvectors.resize(nq);
    }
    coor2_t* p_qvectors = &(qvectors[0].x);
    broadcast(scatter_comm,p_qvectors,3*nq,0);
       
        
    // communicate services to all nodes
    std::string host_str,fileservice_port_str,monitorservice_port_str;
    if (scatter_comm.rank()==0) {
        host_str = boost::asio::ip::host_name();
        fileservice_port_str = boost::lexical_cast<string>(p_hdf5writer->get_endpoint().port());
        monitorservice_port_str = boost::lexical_cast<string>(p_monitorservice->get_endpoint().port());
    }    
    if (world.rank()==0) Info::Inst()->write("Broadcasting service information");	

    broadcast(scatter_comm,host_str,0);
    broadcast(scatter_comm,fileservice_port_str,0);
    broadcast(scatter_comm,monitorservice_port_str,0);
    if (world.rank()==0) {
        Info::Inst()->write(string("Server host name    : ")+host_str);	
        Info::Inst()->write(string("FileService Port    : ")+fileservice_port_str);	
        Info::Inst()->write(string("MonitorService Port : ")+monitorservice_port_str);	
    }


    boost::asio::ip::tcp::resolver resolver(io_service);
    //fileserver
    boost::asio::ip::tcp::resolver::query query(boost::asio::ip::tcp::v4(),host_str,fileservice_port_str,boost::asio::ip::resolver_query_base::numeric_service);
    boost::asio::ip::tcp::resolver::iterator it = resolver.resolve(query);
    boost::asio::ip::tcp::endpoint fileservice_endpoint(*it);
    //monitor
    boost::asio::ip::tcp::resolver::query query2(boost::asio::ip::tcp::v4(),host_str,monitorservice_port_str,boost::asio::ip::resolver_query_base::numeric_service);
    it = resolver.resolve(query2);
    boost::asio::ip::tcp::endpoint monitorservice_endpoint(*it);
    
    if (world.rank()==0) Info::Inst()->write("Setting up parallel environment...");	
    
    // now create scattering device w/ reference to file server
    IScatterDevice* p_ScatterDevice =NULL;
    p_ScatterDevice = ScatterDeviceFactory::create(
    		scatter_comm,
    		sample,
    	    fileservice_endpoint,
            monitorservice_endpoint,
            qvectors);

	world.barrier();	

    // start computation tasks
    if (world.rank()==0) Info::Inst()->write("Starting scattering...");

    if (world.rank()==0) {
        Info::Inst()->write("Resetting progress timer (timer includes stage time)...");   
    }
    if (world.rank()==0) p_monitorservice->timer_reset();
    if (p_ScatterDevice!=NULL) p_ScatterDevice->run();

    world.barrier();
    
    if (world.rank()==0) Info::Inst()->write("Scattering finished...");
    
    // shutdown services
    if (scatter_comm.rank()==0) {
        p_monitorservice->hangup();
        p_hdf5writer->hangup();
    }
    
    world.barrier();

    if (world.rank()==0) print_analysis();
    
	if (world.rank()==0) {
		Info::Inst()->write(string("Aggregating timing information for performance analysis..."));		
	}
	
	timer.stop("total");
	
    Timer performance_timer;
    if (p_ScatterDevice!=NULL) performance_timer = p_ScatterDevice->timer;
    PerformanceAnalyzer perfanal(world,performance_timer); // collect timing information from everybody.
    if (p_ScatterDevice!=NULL) delete p_ScatterDevice;
        
    if (world.rank()==0) {
		perfanal.report();
		perfanal.report_relative(timer.sum("total")*world.size());
    	Info::Inst()->write(string("Total runtime (s): ")+boost::lexical_cast<string>(timer.sum("total")));
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
