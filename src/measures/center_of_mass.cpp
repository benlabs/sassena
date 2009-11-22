/*
 *  center_of_mass.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */
// direct header
#include "measures/center_of_mass.hpp"

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "coor3d.hpp"
#include "sample/coordinate_set.hpp"
#include "control.hpp"

using namespace std;

CenterOfMass::CenterOfMass(Atoms& atoms,Atomselection& cofm_selection,Atomselection& cs_selection, CoordinateSet& cs) {
	
	if (cs.get_representation()!=CARTESIAN) {
        Err::Inst()->write("Center of Mass only implementated for cartesian Coordinate Sets");
        throw;
	}
	
	if (cofm_selection.empty()) {
		Warn::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
		Warn::Inst()->write("Setting Center of mass to (0,0,0)");
		m_center = CartesianCoor3D(0,0,0);
	}

	double xt,yt,zt,m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	size_t iter = 0; 
	size_t cs_iter = 0;
	for(size_t i = 0; i < cofm_selection.booleanarray.size(); ++i)
	{
		if (cofm_selection.booleanarray[i]) {
			if (! cs_selection.booleanarray[i]) {
				Err::Inst()->write("Called CenterOfMass with two incompatible atom selections");
			} 
			mi = atoms[cs_selection[cs_iter]].mass;
			m += mi;
			xt += cs.c1[cs_iter]*mi;
			yt += cs.c2[cs_iter]*mi;
			zt += cs.c3[cs_iter]*mi;
			
			iter++;
		}
		if (cs_selection.booleanarray[i]) cs_iter++;
		
	}
		
	m_center = CartesianCoor3D(xt/m,yt/m,zt/m);
	
}


CenterOfMass::CenterOfMass(Atoms& atoms,Frame& frame,Atomselection& selection) {
	
	if (selection.empty()) {
		cerr << "Warning! Computing Center of Mass for an empty atomselection" << endl;
		cerr << "Setting Center of mass to (0,0,0)" << endl;		
		m_center = CartesianCoor3D(0,0,0);
        return;
	}

	double xt,yt,zt,m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	for (size_t i=0;i<selection.size();i++) {
		mi = atoms[selection[i]].mass;
		m += mi;
		xt += frame.x[selection[i]]*mi;
		yt += frame.y[selection[i]]*mi;
		zt += frame.z[selection[i]]*mi;
	}
	
	m_center =  CartesianCoor3D(xt/m,yt/m,zt/m);
}

CenterOfMass::operator CartesianCoor3D() {
	return m_center;
}

// end of file