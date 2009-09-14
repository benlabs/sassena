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
#include "coordinate_set.hpp"
#include "frame.hpp"
#include "frames.hpp"

// this is our "container", the system, the 'sample'. It contains all the information about the atoms and the time information in form of frames
class Sample {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & atoms;
		ar & frames;
		ar & background;
    }
	/////////////////// 
	
public:	

	Atoms atoms;
	Frames frames;

	// cache values:
	// used for scattering
	double background;

	Sample() { }
	// the sample can be initialized with a system information file: e.g. a pdb
	Sample(std::string filename) : atoms(filename)  { }
	Sample(std::string filename,std::string fileformat) : atoms(filename,fileformat)  { }

	// default routine for reading structure information from file.
	void add_atoms(std::string filename,std::string fileformat) { return atoms.add(filename,fileformat); }
	
	void deuter(std::string group);
	void deuter(std::string group,double coverage,int seed);
	
};

#endif
