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

#ifndef ATOMS_HPP_
#define ATOMS_HPP_

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
#include "atom.hpp"
#include "frame.hpp"

class Atom;
class Frame;

class Atoms : public std::vector<Atom> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<Atom> >(*this);
    }
	///////////////////
public:	
	enum ATOMSFORMAT { PDBFORMAT };

	Atoms() {}
	Atoms(std::string filename, std::string fileformat = "pdb");
	
	void read_particle( std::ifstream& input, std::string fileformat = "pdb");
	void read_solvent( std::ifstream& input, std::string fileformat = "pdb");
	void read_deuter( std::ifstream& input, std::string fileformat = "pdb");

	void write( std::string filename,Frame& frame, std::string fileformat = "pdb");
	void print_statistics();

	void add(std::string filename, std::string fileformat = "pdb");

	~Atoms() {}
	
private:
	std::string guess_atomname(const std::string pdbatomname,std::string fileformat,std::map<std::string,std::string>& quicklookup);

};

#endif
