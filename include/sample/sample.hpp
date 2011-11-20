/** \file
This file contains a the main sample class which defines structure, coordinates and selections.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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

/** 
Managment class for structure and coordinates
*/
class Sample {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & atoms;
		ar & coordinate_sets;

		// have to be set, since coordinate_sets is ignorant
		coordinate_sets.set_atoms(atoms);
		coordinate_sets.set_selection(atoms.selections["system"]);
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
	    
};

#endif

// end of file
