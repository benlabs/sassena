/*
 *  scatter_devices/all_multipole.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef MEASURES__CENTER_OF_MASS_HPP_
#define MEASURES__CENTER_OF_MASS_HPP_

// common header
#include "common.hpp"

// other headers
#include "atoms.hpp"
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "coordinate_set.hpp"


// this helper class takes Atoms and a Coordinateset and returns a cartesian coordinate
class CenterOfMass {
	CartesianCoor3D m_center;
public:
	CenterOfMass(Atoms& atoms,Atomselection& cofm_selection,Atomselection& cs_selection, CoordinateSet& atoms);
	
	operator CartesianCoor3D (); //conversion operator
};


#endif 

// end of file
