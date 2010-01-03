/*
 *  frame.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "sample/frame.hpp"

// standard header

// special library headers

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"

using namespace std;

void Frame::clear() { 
	x.clear(); 
	y.clear(); 
	z.clear(); 
	unitcell.clear();
}

CartesianCoor3D Frame::coord3D(size_t i) { 
	return CartesianCoor3D(x[i],y[i],z[i]);	
}

// end of file

