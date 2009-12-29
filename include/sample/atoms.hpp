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
        ar & system_selection;
    }
	///////////////////
public:	
	
	std::map<std::string,Atomselection> selections;	
    Atomselection system_selection;
    
	Atoms() {}
	Atoms(std::string filename, std::string fileformat = "pdb");
	
	void write( std::string filename,Frame& frame, std::string fileformat = "pdb");

	void add(std::string filename, std::string fileformat = "pdb");

    // support for selections
    std::vector<std::pair<std::string,Atomselection> > select_ndx(std::string filename);
    Atomselection select_pdb(std::string filename, std::string select,double select_value);
    Atomselection select(std::string label="");
    
	void push_selection(std::string name, Atomselection& selection);
    void clear_selections();
	
	void assert_selection(std::string groupname); // throws an exception when groupname isn't found.
	
	~Atoms() {}
};

#endif
