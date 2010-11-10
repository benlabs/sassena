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

#ifndef SAMPLE__CENTER_OF_MASS_HPP_
#define SAMPLE__CENTER_OF_MASS_HPP_

// common header
#include "common.hpp"

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "sample/coordinate_set.hpp"
#include "sample/frame.hpp"

// this helper class takes Atoms and a Coordinateset and returns a cartesian coordinate
class CenterOfMass {
	CartesianCoor3D m_center;
public:
	CenterOfMass(Atoms& atoms,IAtomselection* pcofm_selection,IAtomselection* pcs_selection, CoordinateSet& cs);
	CenterOfMass(Atoms& atoms,Frame& frame,IAtomselection* pselection);
	
	operator CartesianCoor3D (); //conversion operator
};


#endif 

// end of file
