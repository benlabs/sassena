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

void CoordinateSet::rotate(CartesianCoor3D axis1,CartesianCoor3D axis2) {

    // the axis indicates the direction of the new x-axis!
    // we have to rotate the whole system, so that the old x-axis points towards the new one

    CartesianCoor3D al1 = CartesianCoor3D(1,0,0);
    CartesianCoor3D al2 = CartesianCoor3D(0,1,0);
    
    CartesianCoor3D norm1 = axis1/axis1.length();
    CartesianCoor3D norm2 = axis2/axis2.length();

    CartesianCoor3D rot1 = norm1.cross_product(al1); // is : 0 -c b
    double ang1 = acos( norm1*al1 / (norm1.length()*al1.length()) );
    CartesianCoor3D rot2 = norm2.cross_product(al2); // is : 0 -c b
    double ang2 = acos( norm2*al2 / (norm2.length()*al2.length()) );
    

	for(size_t i = 0; i < m_size; ++i)
	{
        double c,s,r_00,r_01,r_02,r_10,r_11,r_12,r_20,r_21,r_22,xo,yo,zo;
        
	    // first alignment
        c = cos(ang1);
        s = sin(ang1);
        r_00 = powf(rot1.x,2) + (1-powf(rot1.x,2))*c;
        r_01 = rot1.x*rot1.y*(1-c) - rot1.z*s;
        r_01 = rot1.x*rot1.z*(1-c) + rot1.y*s;
        r_10 = rot1.x*rot1.y*(1-c) + rot1.z*s;
        r_11 = powf(rot1.y,2) + (1-powf(rot1.y,2))*c;
        r_11 = rot1.y*rot1.z*(1-c) - rot1.x*s;
        r_20 = rot1.x*rot1.z*(1-c) - rot1.y*s;
        r_21 = rot1.y*rot1.z*(1-c) + rot1.x*s;
        r_21 = powf(rot1.z,2) + (1-powf(rot1.z,2))*c;
	    
        xo = x[i];
        yo = y[i];
        zo = z[i];
        	    
		x[i] = r_00*xo + r_01*yo + r_02*zo;
		y[i] = r_10*xo + r_11*yo + r_12*zo;
		z[i] = r_20*xo + r_21*yo + r_22*zo;

        // second alignment
        c = cos(ang2);
        s = sin(ang2);
        r_00 = powf(rot2.x,2) + (1-powf(rot2.x,2))*c;
        r_01 = rot2.x*rot2.y*(1-c) - rot2.z*s;
        r_01 = rot2.x*rot2.z*(1-c) + rot2.y*s;
        r_10 = rot2.x*rot2.y*(1-c) + rot2.z*s;
        r_11 = powf(rot2.y,2) + (1-powf(rot2.y,2))*c;
        r_11 = rot2.y*rot2.z*(1-c) - rot2.x*s;
        r_20 = rot2.x*rot2.z*(1-c) - rot2.y*s;
        r_21 = rot2.y*rot2.z*(1-c) + rot2.x*s;
        r_21 = powf(rot2.z,2) + (1-powf(rot2.z,2))*c;
	    
        xo = x[i];
        yo = y[i];
        zo = z[i];
        	    
		x[i] = r_00*xo + r_01*yo + r_02*zo;
		y[i] = r_10*xo + r_11*yo + r_12*zo;
		z[i] = r_20*xo + r_21*yo + r_22*zo;
		
		
	}	
	
}





// end of file
