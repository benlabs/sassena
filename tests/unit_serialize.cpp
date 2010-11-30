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
#include <boost/date_time.hpp>
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/thread.hpp>

// other headers
#include "control.hpp"
#include "log.hpp"
#include "sample/sample.hpp"

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

    // The rank 0 node is responsible for the progress output and to inform the user
    // compute nodes should be silent all the time, except when errors occur.
    // In that case the size and the rank should be included into the error message in the following way:

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(boost::lexical_cast<string>(string(".Info>>")));
	Warn::Inst()->set_prefix(boost::lexical_cast<string>(string(".Warn>>")));
	Err::Inst()->set_prefix(boost::lexical_cast<string>(string(".Err>>")));
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();

	Sample sample;
	
    bool initstatus = true;
		
    po::variables_map vm;    
    print_title();
    print_description();
    initstatus = init_commandline(argc,argv,vm);


print_initialization();
	    
params->init(vm["config"].as<string>());
database->init(vm["database"].as<string>());	    
sample.init();

detect_parameters();    
read_parameters(vm);
		
	//------------------------------------------//
	//
	// Communication of the sample
	// At this point it is ILLEGAL to change anything within the sample.
	//
	//------------------------------------------//

Info::Inst()->write("Exchanging sample, database & params information with compute nodes... ");
Info::Inst()->write("params... ");
Info::Inst()->write("database... ");
Info::Inst()->write("sample... ");




// params


// database

std::stringstream paramsstream;
{
    boost::archive::text_oarchive ar(paramsstream);
    ar << *params;
}
{
    boost::archive::text_iarchive ar(paramsstream);
    ar >> *params;    
}

std::stringstream databasestream;
{
    boost::archive::text_oarchive ar(databasestream);
    ar << *database;
}
{
    boost::archive::text_iarchive ar(databasestream);
    ar >> *database;    
}

std::stringstream samplestream;
{
    boost::archive::text_oarchive ar(samplestream);
    ar << sample;
}
{
    boost::archive::text_iarchive ar(samplestream);
    ar >> sample;    
}	

Info::Inst()->write("done");

}