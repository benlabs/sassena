/*
 *  atoms.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SAMPLE__ATOMS_HPP_
#define SAMPLE__ATOMS_HPP_

// common header
#include "common.hpp"

// standard header
#include <fstream>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "sample/atom.hpp"
#include "sample/atomselection.hpp"
#include "sample/frame.hpp"

class Atom;
class Frame;
class Atomselection;

class Atoms : public std::vector<Atom> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<Atom> >(*this);
		ar & selections;
    }
	///////////////////
public:	
	
	std::map<std::string,Atomselection> selections;	

	Atoms() {}
	Atoms(std::string filename, std::string fileformat = "pdb");
	
	void read_particle( std::ifstream& input, std::string fileformat = "pdb");
	void read_solvent( std::ifstream& input, std::string fileformat = "pdb");
	void read_deuter( std::ifstream& input, std::string fileformat = "pdb");

	void write( std::string filename,Frame& frame, std::string fileformat = "pdb");
	void print_statistics();

	void add(std::string filename, std::string fileformat = "pdb");
	void add(std::string label);
	
	void add_selection(std::string name, std::string filename, std::string format,std::string select,double select_value);
	void add_selection(std::string name,bool select);
	
	void assert_selection(std::string groupname); // throws an exception when groupname isn't found.
	~Atoms() {}
};

#endif
