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

CoordinateSet::CoordinateSet() {
	m_size = 0;
}

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

CoordinateSet::CoordinateSet(CoordinateSet& original_cs,Atomselection& original_selection, Atomselection& sub_selection) {

	// to construct a new coordinateset from a prior one, the following conditions has to be met:
	// the selection of the new one must be a subset of the old one
	// otherwise an exception is raised. 
	if ( ! sub_selection.is_subset_of(original_selection) ) {
		Err::Inst()->write("Tried to construct a coordinate set from one which it isn't a subset of!");
		throw;
	}
	
	m_size = sub_selection.size();
	x.resize(m_size);
	y.resize(m_size);
	z.resize(m_size);

	// use the complete booleanarray of the selection, this allows translation of the indexes on the fly
	size_t cs_index_iter = 0;
	size_t index_iter = 0;	
	for(size_t i = 0; i < sub_selection.booleanarray.size(); ++i)
	{
		if (sub_selection.booleanarray[i]) {
			x[index_iter] = original_cs.x[cs_index_iter];
			y[index_iter] = original_cs.y[cs_index_iter];
			z[index_iter] = original_cs.z[cs_index_iter];
			index_iter++;
		}		

		if (original_selection.booleanarray[i]) cs_index_iter++;
	}
}

void CoordinateSet::translate(CartesianCoor3D trans,Atomselection& original_selection, Atomselection& sub_selection) {
	// to translate only a subpart, the following conditions has to be met:
	// the sub selection must be a subset of the containing one
	if ( ! sub_selection.is_subset_of(original_selection) ) {
		Err::Inst()->write("Tried to translate a coordinate with an incompatible atomselection!");
		throw;
	}

	size_t index_iter = 0;	
	for(size_t i = 0; i < original_selection.booleanarray.size(); ++i)
	{
		if (sub_selection.booleanarray[i]) {
			x[index_iter] += trans.x;
			y[index_iter] += trans.y;
			z[index_iter] += trans.z;
		}		

		if (original_selection.booleanarray[i]) index_iter++;
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
