/*
 *  coordinateset.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef COORDINATESET_HPP_
#define COORDINATESET_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "frame.hpp"

//forward declaration...

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

// This class is used by Frame to store selection specific coordinate arrays
class CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & x;
		ar & y;
		ar & z;
		ar & m_size;
    }

	size_t m_size;
public:
	
	CoordinateSet(Frame& frame,Atomselection& selection); 
	
	std::vector<double> x; // x-coordinates
	std::vector<double> y; // y-coordinates
	std::vector<double> z; // z-coordinates
		
	size_t size() { return m_size; }

	void translate(CartesianCoor3D trans);
};

#endif

// end of file