/*
 *  sample.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "sample.hpp"

// standard header

// special library headers

// other headers


using namespace std;

void Sample::deuter(std::string group) {
	Atomselection as = atomselections[group];
	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		atoms[*asi].name="deuterium";		
	}
}

void Sample::add_frame(std::string filename,string filetype) { 

	if (filetype=="dcd") {
		dcdframes.add_file(filename,atoms);		
	}
}

// end of file
