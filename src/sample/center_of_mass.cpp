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
#include "sample/center_of_mass.hpp"

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "sample/coordinate_set.hpp"
#include "control.hpp"
#include "log.hpp"

using namespace std;

CenterOfMass::CenterOfMass(Atoms& atoms,IAtomselection* pcofm_selection,IAtomselection* pcs_selection, CoordinateSet& cs) {
	
	IAtomselection& cs_selection = *pcs_selection;
    IAtomselection& cofm_selection = *pcofm_selection;
    
	
	if (cs.get_representation()!=CARTESIAN) {
        Err::Inst()->write("Center of Mass only implementated for cartesian Coordinate Sets");
        throw;
	}
	
    size_t csel_total = cs_selection.size();
    size_t ssel_total = cofm_selection.size();
    
	if (ssel_total<1) {
		Warn::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
		Warn::Inst()->write("Setting Center of mass to (0,0,0)");
		m_center = CartesianCoor3D(0,0,0);
	} else if (csel_total<1) {
		Warn::Inst()->write("Warning! Computing Center of Mass for an empty coordinate set");
		Warn::Inst()->write("Setting Center of mass to (0,0,0)");
		m_center = CartesianCoor3D(0,0,0);	    
	} else {

        coor2_t xt,yt,zt;
        double m,mi;
    	xt = yt = zt = 0.0;
    	m = 0.0;

        size_t csel_iter = 0;
        size_t ssel_iter = 0;

        while( ( csel_iter < csel_total) && (ssel_iter < ssel_total) ) {
            size_t csel_index = cs_selection[csel_iter];
            size_t ssel_index = cofm_selection[ssel_iter];

            if (csel_index==ssel_index) {
                size_t csel_id = atoms[csel_index];
                mi = Database::Inst()->masses.get(csel_id);
    			m += mi;
    			
    			xt += cs.c1[csel_iter]*mi;
    			yt += cs.c2[csel_iter]*mi;
    			zt += cs.c3[csel_iter]*mi;
    			
                ssel_iter++;
                csel_iter++;           
            } else if (csel_index>ssel_index) {
                ssel_iter++;
            } else if (csel_index<ssel_index) {
                csel_iter++;
            }
        }
        
    	m_center = CartesianCoor3D(xt/m,yt/m,zt/m);

	}		
	
}


CenterOfMass::CenterOfMass(Atoms& atoms,Frame& frame,IAtomselection* pselection) {
	
    IAtomselection& selection = *pselection;
	
	if (selection.size()==0) {
		Err::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
		Err::Inst()->write("Setting Center of mass to (0,0,0)");		
		m_center = CartesianCoor3D(0,0,0);
        return;
	}

    coor2_t xt,yt,zt;
    double m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	for (size_t i=0;i<selection.size();i++) {
        size_t sel_id = atoms[selection[i]];
        mi = Database::Inst()->masses.get(sel_id);
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
