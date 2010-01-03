/*
 *  particletrajectory.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef PARTICLETRAJECTORY_HPP_
#define PARTICLETRAJECTORY_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <fstream>
#include <string>
#include <vector>

// special library headers
#include <boost/regex.hpp> 
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/frame.hpp"

class ParticleTrajectory {
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & x;
		ar & y;
		ar & z;
		
		ar & m_size;
		ar & m_atom_selection_index;
    }
	///////////////////	

	Atomselection* p_selection;

	size_t m_size;	
	size_t m_atom_selection_index; // corresponds to the index position in the selection
public: 	
	std::vector<double> x,y,z;

	ParticleTrajectory() : m_size(0), m_atom_selection_index(0) { }
	ParticleTrajectory(size_t asi) : m_size(0), m_atom_selection_index(asi) {}
	
	CartesianCoor3D coor(size_t i) { return CartesianCoor3D(x[i],y[i],z[i]);}
	
	void push_back(CartesianCoor3D r) { x.push_back(r.x); y.push_back(r.y); z.push_back(r.z); m_size++;}
	size_t atom_selection_index() { return m_atom_selection_index; }
	size_t size() { return m_size;}
};

#endif

// end of file
