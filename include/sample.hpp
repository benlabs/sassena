/*
 *  sample.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SAMPLE_HPP_
#define SAMPLE_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <string>

// special library headers
#include <boost/serialization/access.hpp>

// other headers
#include "atoms.hpp"
#include "atomselection.hpp"
#include "atomselections.hpp"
#include "frame.hpp"
#include "frames.hpp"

// this is our "container", the system, the 'sample'. It contains all the information about the atoms and the time information in form of frames
class Sample {
	
public:	
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & atoms;
        ar & atomselections;
        ar & dcdframes;
		ar & framecache;
		ar & curframe_i;
		ar & framecache_max;
		ar & wrapping;
		ar & centergroup;
		ar & background;
    }
	/////////////////// 
	
	Atoms atoms;
	Atomselections atomselections;
	DcdFrames dcdframes;
		
	std::map<int,DcdFrame> framecache;
	
	int framecache_max;
	int curframe_i;
	// the currentframe is a cache, a buffer 
	// this enables the analysis routines to work on a default target

	// unit cell behaviour:
	bool wrapping;
	std::string centergroup;

	// cache values:
	// used for scattering
	double background;

	Sample() { framecache_max = 2;}
	// the sample can be initialized with a system information file: e.g. a pdb
	Sample(std::string filename) : atoms(filename)  { framecache_max = 2;}
	Sample(std::string filename,std::string fileformat) : atoms(filename,fileformat)  { framecache_max = 2; }

	void add_frame(std::string filename,std::string filetype);

	void add_selection(std::string name, std::string filename, std::string format,std::string select,double select_value) {
		atomselections[name] = Atomselection(atoms,filename,format,select,select_value);
	}
	void add_selection(std::string name, std::string filename, std::string format) {
		atomselections[name] = Atomselection(atoms,filename,format);		
	}
	void add_selection(std::string name, Atoms::iterator b, Atoms::iterator e) {
		atomselections[name] = Atomselection(b,e);
	}
	
	// default routine for reading structure information from file.
	void add_atoms(std::string filename) { return atoms.add(filename); }
	
	void read_frame(int framenumber) { 
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
		}
		curframe_i = framenumber;
	}
	
	DcdFrame& currentframe() { 
		return framecache[curframe_i];	
	}
	
	void deuter(std::string group);
	
	void write_frame(std::string filename, std::string af = "pdb") { atoms.write(filename,currentframe(),af); }	
};

#endif
