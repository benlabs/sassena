/*
 *  frame.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SAMPLE__FRAME_HPP_
#define SAMPLE__FRAME_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"

//forward declaration...
class Atom;
class Atoms;

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

class Frame {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & number_of_atoms;
		ar & t;
		ar & origin;
		ar & x;
		ar & y;
		ar & z;
		ar & unitcell;
    }
	
public:
	size_t number_of_atoms; //store this for performance
	
	double t; // time information
	
	CartesianCoor3D origin; // this is
	
	std::vector<coor2_t> x; // x-coordinates
	std::vector<coor2_t> y; // y-coordinates
	std::vector<coor2_t> z; // z-coordinates
	
	// base vectors of the unit cell
	std::vector<CartesianCoor3D> unitcell;
	
	void clear(); // make frame 'empty'

	CartesianCoor3D coord3D(size_t i); // get coordinates for i'th atom
};

#endif

// end of file
