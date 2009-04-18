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
#include "settings.hpp"

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

void Sample::read_frame(int framenumber) { 
	if (framecache.find(framenumber)!=framecache.end()) {
		curframe_i = framenumber;
	}
	else {
		if (framecache.size()>framecache_max) {
			framecache.clear();
		}
		DcdFrame& cf = framecache[framenumber];
		cf.clear();
		dcdframes.read(framenumber,cf); 
		if (wrapping) {
			cf.origin = cf.cofm(atoms,atomselections[centergroup]);
			cf.wrap(); 					
		}	
		string target = Settings::get("main")["scattering"]["target"];
		Atomselection& ast = atomselections[target];
		dcdframes.extracttarget(cf,ast); // generates a 3xN matrix within the frame class to hold only target atoms	
	}
	curframe_i = framenumber;
}

// end of file
