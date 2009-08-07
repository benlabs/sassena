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
#include "log.hpp"
#include "parameters.hpp"
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
	
//	boost::mpi::environment env(argc, argv);
// 	boost::mpi::communicator world;
//	int rank = world.rank();
//	int size = world.size();
  	
	
	clog << "INFO>> " << "Reading configuration file " << argv[1] << endl;
	// read in configuration file, delete command line arguments
	if (!Settings::read(argc,argv)) { cerr << "Error reading the configuration" << endl; throw; }

	// read structure
	Sample sample;
	sample.add_atoms(Settings::get_filepath(Settings::get("main")["sample"]["structure"]["file"]), Settings::get("main")["sample"]["structure"]["format"]);
	
	// read frames
//	for (int i=0;i<Settings::get("main")["sample"]["frames"].getLength();i++) {
//		clog << "INFO>> " << "Reading frames from: " << (const char *) Settings::get("main")["sample"]["frames"][i]["file"] << endl;
 //       sample.frames.add_frameset(Settings::get("main")["sample"]["frames"][i]["file"],Settings::get("main")["sample"]["frames"][i]["type"],sample.atoms);
//	}
			
	// estimate rank logic 
	
	std::vector<CartesianCoor3D> qqqvectors;
	std::vector<int> frames;
	Settings::get_qqqvectors(qqqvectors);
	// fill in frames
	int framestride = 1;
	if (Settings::get("main")["scattering"].exists("framestride")) {
		framestride = Settings::get("main")["scattering"]["framestride"];
	}
	for(size_t i=0;i<sample.frames.size();i+= framestride) {
		frames.push_back(i);
	}	
	
	size_t nf = frames.size();
	size_t nq = qqqvectors.size();

	clog << "INFO>> Best choice for NN(number of nodes): " << nf << endl;
	for(size_t i = 1; i < nf; ++i)
	{
		if ((nf % i)==0) {
			clog << "INFO>> Overloading choice: ";		
			clog << i << " ("<< nf / i << "ZE)" << endl;
		}
	}
		
	
	return 0;
}

// end of file
