/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
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
	CenterOfMass(Atoms& atoms, CoordinateSet& cs, IAtomselection* pcs_selection, IAtomselection* pcofm_selection);
	CenterOfMass(Atoms& atoms,Frame& frame,IAtomselection* pselection);
	CenterOfMass(Atoms& atoms,CoordinateSet& cs,IAtomselection* pselection);
	
	operator CartesianCoor3D (); //conversion operator
};

// operational class, modifies arguments
class Fit {
public:
    Fit(
        Atoms& atoms, // contains IDs in sequence
        CoordinateSet& cs, // contains original coordinates of target
        IAtomselection* pcs_selection, // selection corresponding to target
        IAtomselection* pcs_selection_manip, // (sub-)selection of target which specifies which atoms to move
        CoordinateSet& cs_ref,  // coordinate set of the reference structure
        IAtomselection* pref_selection ); // selection correspoding to reference    
};


#endif 

// end of file
