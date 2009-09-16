/*
 *  coordinateset.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */
// direct header
#include "coordinate_set.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "frame.hpp"
#include "coor3d.hpp"
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"

using namespace std;

CoordinateSet::CoordinateSet(Frame& frame,Atomselection& selection) {

	m_size = selection.size();
	x.resize(m_size);
	y.resize(m_size);
	z.resize(m_size);

	for(size_t i = 0; i < m_size; ++i)
	{
		size_t thisindex = selection[i] ;
		x[i] = frame.x[thisindex];
		y[i] = frame.y[thisindex];
		z[i] = frame.z[thisindex];		
	}
}

void CoordinateSet::translate(CartesianCoor3D trans) {

	for(size_t i = 0; i < m_size; ++i)
	{
		x[i] += trans.x;
		y[i] += trans.y;
		z[i] += trans.z;	
	}	
	
}

// end of file
