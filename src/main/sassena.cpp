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
#include "report/progress_reporter.hpp"
#include "report/timer.hpp"
#include "sample/sample.hpp"
#include "scatter_devices/scatter_device_factory.hpp"
#include "scatter_devices/scatter_spectrum.hpp"
#include "threadpool/threadpool.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;
bool init_flag = false;
boost::asio::ip::tcp::endpoint server_endpoint;

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

        ("scattering-data-file",po::value<string>()->default_value("fqt.h5"),"name of the data output file")
        ("limits-computation-threads",po::value<string>()->default_value("1"),"threadpool size for orientational averaging")        
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    
    if (vm.find("help")!=vm.end()) {
        cout << desc << endl;
        return false;
    }
    
    if (vm.find("config")==vm.end()) {
        Info::Inst()->write("Require configuration file");
        cout << desc << endl;
        return false;
    }
        
    if (vm.find("database")==vm.end()) {
        Info::Inst()->write("Require database file");
        cout << desc << endl;
        return false;
    }
    
    return true;
}

void read_parameters(po::variables_map vm) {
    Info::Inst()->write("Checking command line for parameter overwrite");
    // second stage , allow command line parameters to overwrite defaults in the configuration file
	if (!vm["limits-computation-threads"].defaulted()) {
        Info::Inst()->write(string("Overwrite: limits.computation.threads=")+vm["limits-computation-threads"].as<string>());
        Params::Inst()->limits.computation.threads = boost::lexical_cast<size_t>(vm["limits-computation-threads"].as<string>());
    }
    if (!vm["scattering-data-file"].defaulted()) {
        Info::Inst()->write(string("Overwrite: scattering.data.file=")+vm["scattering-data-file"].as<string>());
        Params::Inst()->scattering.data.file = vm["scattering-data-file"].as<string>();        
    }
}

void detect_parameters() {
    if (boost::thread::hardware_concurrency()!=0) {
        Info::Inst()->write(string("Detect: Number of Processors per machine: ")+boost::lexical_cast<string>(boost::thread::hardware_concurrency()));   
    }
}

void monitor_server_thread(size_t NN) {
    // setup monitoring service
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket( io_service );
    boost::asio::ip::tcp::acceptor acceptor(io_service,server_endpoint);

    // write port into the reference variable
    server_endpoint = acceptor.local_endpoint();
    // notify successful initialization
    init_flag = true;
    
    stringstream ep_strstr; ep_strstr << server_endpoint; 
    Info::Inst()->write(string("S: Monitoring service setup at ")+ep_strstr.str());

    vector<double> progresses(NN);
    
    ProgressReporter progress_reporter(0,NN); // this acts like a progress bar
    double maxprogress = NN*1.0;
    double totalprogress=0;
    
    while (totalprogress!=maxprogress) {
        acceptor.accept(socket);
        
        size_t rank; 
        double progress;

        socket.read_some(boost::asio::buffer(&rank,sizeof(size_t)));
        socket.read_some(boost::asio::buffer(&progress,sizeof(double)));

        progresses[rank]=progress;
        
        double sum=0;
        for(size_t i=0;i<NN;i++) {
            sum+= progresses[i];
        }
        totalprogress = sum;
        progress_reporter.set(totalprogress);
		progress_reporter.report(); // only reports if significant progress has been made!
		
        socket.close();
    }
}

void monitor_client_thread(boost::mpi::communicator scatter_comm,ScatterDevice* pScatterDevice,string host_str,string port_str) {
    // setup monitoring service
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket( io_service );
    
    boost::asio::ip::tcp::resolver resolver(io_service);
    boost::asio::ip::tcp::resolver::query query(boost::asio::ip::tcp::v4(),host_str,port_str);
    boost::asio::ip::tcp::resolver::iterator it = resolver.resolve( query );

    boost::asio::ip::tcp::endpoint endpoint(*it);

    size_t rank = scatter_comm.rank();
    double progress = 0.0;
    while (progress!=1.0) {

        if (progress!=pScatterDevice->progress()) {
            progress = pScatterDevice->progress();
            stringstream ep_strstr2; 
            ep_strstr2 << endpoint; 
            Info::Inst()->write(string("C: Monitoring service setup at ")+ep_strstr2.str());
            
            socket.connect(*it);
            socket.write_some(boost::asio::buffer(&rank,sizeof(size_t)));
            socket.write_some(boost::asio::buffer(&progress,sizeof(double)));
            socket.close();
        }
		
		sleep(2);        
    }
}

void scatter_thread(boost::mpi::communicator scatter_comm,ScatterDevice* pScatterDevice) {
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        if (scatter_comm.rank()==0) {		
    		Info::Inst()->write("Qvectors orientations used for averaging: ");
            for (size_t i=0;i<Params::Inst()->scattering.average.orientation.vectors.size();i++) {
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
	
    if (scatter_comm.rank()==0) {		
		Info::Inst()->write("Starting scattering calculation");
	}
	
    while(pScatterDevice->status()==0) {
		pScatterDevice->compute();
		pScatterDevice->write();
        pScatterDevice->next();
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
	
	Info::Inst()->set_prefix(to_s(world.rank())+string(".Info>>"));
	Warn::Inst()->set_prefix(to_s(world.rank())+string(".Warn>>"));
	Err::Inst()->set_prefix(to_s(world.rank())+string(".Err>>"));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
	Timer timer;
	timer.start("total");

    bool initstatus = true;
		
    po::variables_map vm;    
    if (world.rank()==0) {
        print_title();
        initstatus = init_commandline(argc,argv,vm);
    }
    
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
    
	if (world.rank()==0) Info::Inst()->write(string("Set background scattering length density set to ")+to_s(Params::Inst()->scattering.background.factor));
		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

	if (world.rank()==0) Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
    
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

    if (world.rank()==0) print_computation();

    boost::mpi::communicator scatter_comm = boost::mpi::communicator(world,world.group());

	world.barrier();

    if (world.rank()==0) Info::Inst()->write("Setting up parallel environment...");	
    	    
    	    
    ScatterDevice* p_ScatterDevice =NULL;
    p_ScatterDevice = ScatterDeviceFactory::create(
    		scatter_comm,
    		sample,
    		params->scattering.qvectors );

	world.barrier();



        // setup monitoring service
    //    boost::asio::io_service io_service;
    //    boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::tcp::v4(),11111);
    //    boost::asio::ip::tcp::acceptor acceptor(io_service,endpoint);
//    stringstream port_strstr; port_strstr <<  acceptor.local_endpoint().port();
        string host_str = boost::asio::ip::host_name();
    //    boost::mpi::broadcast(scatter_comm,port_str,0);
    //    boost::mpi::broadcast(scatter_comm,host_str,0);

    
    
    string port_str = "11111";
    boost::thread* p_pthread = NULL;
    if (world.rank()==0) {
        // create a lock
        boost::mutex init_mutex;
        boost::condition init_condition;

        Info::Inst()->write(string("Setting up monitoring service..."));
        boost::mutex::scoped_lock init_lock( init_mutex );
        boost::asio::ip::tcp::endpoint endpoint(boost::asio::ip::tcp::v4(),0);
        
        p_pthread = new boost::thread(boost::bind(&monitor_server_thread,scatter_comm.size()));        
        stringstream ep_strstr; 

        while (init_flag!=true) {
            if (init_flag) init_lock.unlock();
            stringstream ep_strstr2; 
            ep_strstr2 << server_endpoint; 
            Info::Inst()->write(string("W: Monitoring service setup at ")+ep_strstr2.str());
            Info::Inst()->write(string("W: Init flag= ")+boost::lexical_cast<string>(init_flag));

            init_condition.timed_wait(init_lock,boost::posix_time::microseconds(1000000));            
        }
        ep_strstr << server_endpoint;
        
        Info::Inst()->write(string("Monitoring service setup at ")+ep_strstr.str());
    }

    port_str = boost::lexical_cast<string>(server_endpoint.port());
    if (world.rank()==0) Info::Inst()->write(string("Broadcasting ")+host_str + string(":") + port_str + string(" to clients"));

    boost::mpi::broadcast(world,host_str,0);    
    boost::mpi::broadcast(world,port_str,0);
    
    // start computation tasks
    if (world.rank()==0) Info::Inst()->write("Starting computation...");

    boost::thread pthread(boost::bind(&monitor_client_thread,scatter_comm,p_ScatterDevice,host_str,port_str));    
    boost::thread sthread(boost::bind(&scatter_thread,scatter_comm,p_ScatterDevice));

    // do a wait on the threads here
    sthread.join();
    pthread.timed_join(boost::posix_time::seconds(1));
    pthread.interrupt();
    if (p_pthread!=NULL) {
        p_pthread->timed_join(boost::posix_time::seconds(1));
        p_pthread->interrupt();
        delete p_pthread;
    }

    if (world.rank()==0) print_analysis();
    
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
