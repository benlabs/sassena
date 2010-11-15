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

class Atoms {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar.register_type(static_cast<IndexAtomselection*>(NULL));
        ar.register_type(static_cast<RangeAtomselection*>(NULL));
        
		ar & selections;
        ar & ids_;
    }
	///////////////////
	
    std::vector<size_t> ids_;
public:	
	
	std::map<std::string,IAtomselection*> selections;	
    
	Atoms() {}
	Atoms(std::string filename, std::string fileformat = "pdb");
	
//	void write( std::string filename,Frame& frame, std::string fileformat = "pdb");

	void add(std::string filename, std::string fileformat = "pdb");
    
    IAtomselection* select(std::string expression);
	void push_selection(std::string name, IAtomselection* selection);
    void clear_selections();
	
	void assert_selection(std::string groupname); // throws an exception when groupname isn't found.
	
    ~Atoms();
    
    size_t operator[](size_t index) { return ids_[index]; }
    size_t size() { return ids_.size(); }    
};

#endif

// end of file
