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
#include "control.hpp"
#include "sample.hpp"


#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
using namespace std;

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
		Info::Inst()->write("1. Sassena - Elastic Scattering Calculations on Parallel Computers       ");
		Info::Inst()->write("   to be published                                                       ");		
		Info::Inst()->write(".........................................................................");
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
	
	    
	    Info::Inst()->write("reading params");
        params->init(string(argv[1]));

        Info::Inst()->write("reading database");
        database->init(string(argv[2]));
	    
        sample.init();
	    

		//------------------------------------------//
		//
		// Preparation, Analysis of the system
		//
		//------------------------------------------//
	
		//------------------------------------------//
		//
		// Communication of the sample
		// At this point it is ILLEGAL to change anything within the sample.
		//
		//------------------------------------------//

		// before calculating anything we need to communicate the sample to any node
		// this is a one-to-all communication

        std::string sample_filename = string("dumpcontrol-sample-")+to_s(world.rank())+string(".bin");
        std::string params_filename = string("dumpcontrol-params-")+to_s(world.rank())+string(".bin");
        std::string database_filename = string("dumpcontrol-database-")+to_s(world.rank())+string(".bin");

		std::ofstream sample_ofile(sample_filename.c_str());
		std::ofstream database_ofile(database_filename.c_str());
		std::ofstream params_ofile(params_filename.c_str());

		sample.coordinate_sets.clear_cache(); // reduce overhead

		boost::archive::text_oarchive sample_oarchive(sample_ofile);
		boost::archive::text_oarchive database_oarchive(database_ofile);
		boost::archive::text_oarchive params_oarchive(params_ofile);


		params_oarchive << *params ;
		database_oarchive << *database ;
		sample_oarchive << sample ;

		sample_ofile.close();
		database_ofile.close();
		params_ofile.close();
		

	return 0;
}

// end of file
