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
#include "mpi/wrapper.hpp"
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

    std::stringstream paramsstream;
    if (world.rank()==0) {
        boost::archive::text_oarchive ar(paramsstream); 
        ar << *params;
    }
	mpi::wrapper::broadcast_stream(world,paramsstream,0);
    if (world.rank()!=0) {
        boost::archive::text_iarchive ar(paramsstream); 
        ar >> *params;
    }

	world.barrier();

    if (world.rank()==0) Info::Inst()->write("database... ");

    std::stringstream databasestream;
    if (world.rank()==0) {
        boost::archive::text_oarchive ar(databasestream); 
        ar << *database;
    }
	mpi::wrapper::broadcast_stream(world,databasestream,0);
    if (world.rank()!=0) {
        boost::archive::text_iarchive ar(databasestream); 
        ar >> *database;
    }

	world.barrier();
    
    if (world.rank()==0) Info::Inst()->write("sample... ");

    std::stringstream samplestream;
    if (world.rank()==0) {
        boost::archive::text_oarchive ar(samplestream); 
        ar << sample;
    }
	mpi::wrapper::broadcast_stream(world,samplestream,0);
    if (world.rank()!=0) {
        boost::archive::text_iarchive ar(samplestream); 
        ar >> sample;
    }
    
	world.barrier();
	
	timer.stop("sample::communication");

	if (world.rank()==0) Info::Inst()->write("done");

}