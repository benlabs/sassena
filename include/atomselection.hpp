/*
 *  atomselection.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef ATOMSELECTION_HPP_
#define ATOMSELECTION_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

// other headers
#include "atom.hpp"
#include "atoms.hpp"
class Atoms;

class Atomselection : public std::vector<size_t> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<size_t> >(*this);
		ar & booleanarray;
		ar & name;
    }
	///////////////////
public:	
	std::vector<bool> booleanarray; // map an atom index position to a yes/no value
	std::string name;
	
	enum PDBSELECT {SEGID,RESID,RESNAME,BETA};

	Atomselection(Atoms& atoms,std::string filename,std::string format,std::string select,double select_value,std::string name="");
	Atomselection(Atoms& atoms,std::string filename,std::string format,std::string name="");
	Atomselection(Atoms& atoms,bool select = false,std::string name="");
	Atomselection(Atoms& atoms,std::vector<size_t> indexes,std::string name="");

	Atomselection& operator+=(const Atomselection& );
	// empty atomselection is actually invalid!
	Atomselection() {}
	
	void add(const Atom& atom);
	void remove(const Atom& atom);
};

#endif
