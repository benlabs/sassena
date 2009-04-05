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
#include <list>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

// other headers
#include "atom.hpp"
#include "atoms.hpp"

class Atomselection : public std::list<int> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::list<int> >(*this);
    }
	///////////////////
public:	
	enum PDBSELECT {SEGID,RESID,RESNAME,BETA};

	Atomselection(Atoms& atoms,std::string filename,std::string format,std::string select,double select_value);
	Atomselection(Atoms& atoms,std::string filename,std::string format);
	// selects by index+length
	Atomselection(Atoms::iterator f,int length);
	// selects by range
	Atomselection(Atoms::iterator f, Atoms::iterator l);
	// select my "atomname" match
	Atomselection(Atoms::iterator f,Atoms::iterator l, std::string match);
	// select whole set of atoms, w/ match
//	Atomselection(std::string match);
	// select whole set of atoms	
	Atomselection(Atoms& atoms);
	Atomselection() {};
	
	void add(Atom& atom);
	void readpdb(Atoms& atoms,std::string pdbname,PDBSELECT select);
		
	// 
//	void truncate_sphere(CartesianCoor3D origin,double radius);
//	void truncate_sshell(CartesianCoor3D origin,double radiusinner,double radiusouter);

//	void truncate_cylinder(Sample& sample,CartesianCoor3D origin,CartesianCoor3D zaxis,double radius);
//	void truncate_cshell(Sample& sample,CartesianCoor3D origin,CartesianCoor3D zaxis,double radiusinner,double radiusouter);
	
};

#endif
