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
#include <libconfig.h++>
#include <boost/math/special_functions.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "sample.hpp"
#include "scatterdata.hpp"
#include "settings.hpp"
#include "tasks.hpp"

using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace boost::filesystem;
namespace mpi = boost::mpi;

int main (int argc, char **argv)
{
  	// sassconf is a helper to build configuration files on the fly
	// it provides 3 functionalities:
	// generate
	// set
	// get
	
	// generate [default]
	// is used to construct a default configuration file
	// an additional argument can be used to describe various "defaults" 
	
	// set path type value [value]
	// is used to set a specific configuration parameter to the given value
	// there are some generic values which are protected. These correspond to the built-in datatypes: double,int,string,array/vector,list,map/hash
	
	// get path
	// retrieves values 
	
	clog << "INFO>> " << "Reading configuration file " << argv[1] << endl;
	// read in configuration file, delete command line arguments
	if (!Settings::read(argc,argv)) { cerr << "Error reading the configuration" << endl; throw; }

	// read structure
	Sample sample;
	sample.add_atoms(Settings::get_filepath(Settings::get("main")["sample"]["structure"]["file"]), Settings::get("main")["sample"]["structure"]["format"]);
	
	// read frames
	for (int i=0;i<Settings::get("main")["sample"]["frames"].getLength();i++) {
		clog << "INFO>> " << "Reading frames from: " << (const char *) Settings::get("main")["sample"]["frames"][i]["file"] << endl;
        sample.frames.add_frameset(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"],sample.atoms);
	}

	
	return 0;
}

//end of file
