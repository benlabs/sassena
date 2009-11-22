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

#ifndef SAMPLE__SAMPLE_HPP_
#define SAMPLE__SAMPLE_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <string>

// special library headers
#include <boost/serialization/access.hpp>

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/coordinate_set.hpp"
#include "sample/coordinate_sets.hpp"


// this is our "container", the system, the 'sample'. It contains all the information about the atoms and the time information in form of frames
class Sample {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & atoms;
		ar & coordinate_sets;
    }
	/////////////////// 
	
public:	

	Atoms atoms;
    CoordinateSets coordinate_sets;

    Sample();
	// the sample can be initialized with a system information file: e.g. a pdb
	Sample(std::string filename) : atoms(filename)  { }
	Sample(std::string filename,std::string fileformat) : atoms(filename,fileformat)  { }

    void init();

	// default routine for reading structure information from file.
	void add_atoms(std::string filename,std::string fileformat) { return atoms.add(filename,fileformat); }
	
	void deuter(std::string group);
	void deuter(std::string group,double coverage,int seed);
    
    void set_kappa(std::string group, double kappa);
    
};

#endif

// end of file
