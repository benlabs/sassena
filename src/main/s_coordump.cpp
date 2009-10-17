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
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "decomposition_plan.hpp"
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"
#include "performance_analyzer.hpp"
#include "progress_reporter.hpp"
#include "sample.hpp"
#include "scatter_devices.hpp"
#include "scatter_spectrum.hpp"
#include "timer.hpp"

using namespace std;

int main(int argc,char** argv) {

	Sample sample;

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	Params::Inst();
	
	Info::Inst()->set_prefix(string("Info>>"));
	Warn::Inst()->set_prefix(string("Warn>>"));
	 Err::Inst()->set_prefix(string("Err>>"));
	
	if (argc!=4) {
        Err::Inst()->write("Need 3 arguments: parameters database outputfile");
        throw;
	}
	
	Params* params = Params::Inst();
	Database* database = Database::Inst();
	
	Info::Inst()->write("reading params");
    params->init(string(argv[1]));

    Info::Inst()->write("reading database");
    database->init(string(argv[2]));
	    
    sample.init();

    CoordinateSets csets;
    csets.set_sample(sample);
    csets.set_selection(sample.atoms.selections["system"]);
    
    for(size_t i = 0; i < sample.frames.size(); ++i) {
        csets.load(i);
    }
    
    if ( boost::filesystem::exists(string(argv[3])) ) {
        Err::Inst()->write(string("Output file ")+string(argv[3])+string(" already exists. Please delete manually."));
        throw;
    }
    csets.write_xyz(string(argv[3]));
    
		
	return 0;
}

// end of file
