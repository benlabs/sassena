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
#include "sample/coordinate_set.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "sample/frame.hpp"
#include "math/coor3d.hpp"
#include "control.hpp"


using namespace std;

CoordinateSet::CoordinateSet() {
	m_size = 0;
}


CoordinateSet::CoordinateSet(CoordinateSet& cs,Atomselection& cs_selection, Atomselection& sub_selection) {
    Atomselection& csel = cs_selection;
    Atomselection& ssel = sub_selection;

    size_t csel_total = csel.indexes.size();
    size_t ssel_total = ssel.indexes.size();
    
    if (csel_total<1) {
        return;
    }
    if (ssel_total<1) {
        return;
    }
    
    size_t csel_iter = 0;
    size_t ssel_iter = 0;
    size_t count =0;
    while( ( csel_iter < csel_total) && (ssel_iter < ssel_total) ) {
        size_t csel_index = csel.indexes[csel_iter];
        size_t ssel_index = ssel.indexes[ssel_iter];
        
        if (csel_index==ssel_index) {
    		c1.push_back(cs.c1[csel_iter]);
    		c2.push_back(cs.c2[csel_iter]);
    		c3.push_back(cs.c3[csel_iter]);
            count++;
            
            ssel_iter++;
            csel_iter++;           
        } else if (csel_index>ssel_index) {
            ssel_iter++;
        } else if (csel_index<ssel_index) {
            csel_iter++;
        }
    }
    
    m_size = count;
    m_representation = cs.get_representation();
}

CartesianCoordinateSet::CartesianCoordinateSet(CartesianCoordinateSet& cs,Atomselection& cs_selection, Atomselection& sub_selection) :
 CoordinateSet(cs,cs_selection,sub_selection)
{
    m_representation = CARTESIAN;
    // inherit copy constructor from parent
}

CartesianCoordinateSet::CartesianCoordinateSet(Frame& frame,Atomselection& selection) {
    m_representation = CARTESIAN;
    
    vector<double>& x = c1;
    vector<double>& y = c2;
    vector<double>& z = c3;

	m_size = selection.indexes.size();
	x.resize(m_size);
	y.resize(m_size);
	z.resize(m_size);

	for(size_t i = 0; i < m_size; ++i)
	{
		size_t thisindex = selection.indexes[i];
		x[i] = frame.x[thisindex];
		y[i] = frame.y[thisindex];
		z[i] = frame.z[thisindex];		
	}
}

void CartesianCoordinateSet::translate(CartesianCoor3D trans, Atomselection& cs_selection, Atomselection& sub_selection) {
	
    Atomselection& csel = cs_selection;
    Atomselection& ssel = sub_selection;

    size_t csel_total = csel.indexes.size();
    size_t ssel_total = ssel.indexes.size();
    
    if (csel_total<1) {
        return;
    }
    if (ssel_total<1) {
        return;
    }
    
    size_t csel_iter = 0;
    size_t ssel_iter = 0;
    
    while( ( csel_iter < csel_total) && (ssel_iter < ssel_total) ) {
        size_t csel_index = csel.indexes[csel_iter];
        size_t ssel_index = ssel.indexes[ssel_iter];
        
        if (csel_index==ssel_index) {
			c1[csel_iter] += trans.x;
			c2[csel_iter] += trans.y;
			c3[csel_iter] += trans.z;
			
            ssel_iter++;
            csel_iter++;           
        } else if (csel_index>ssel_index) {
            ssel_iter++;
        } else if (csel_index<ssel_index) {
            csel_iter++;
        }
    }
    
}


void CartesianCoordinateSet::translate(CartesianCoor3D trans) {
	for(size_t i = 0; i < m_size; ++i)
	{
		c1[i] += trans.x;
		c2[i] += trans.y;
		c3[i] += trans.z;	
	}	
	
}

void CartesianCoordinateSet::rotate(CartesianCoor3D axis1,CartesianCoor3D axis2) {

    // the axis indicates the direction of the new x-axis!
    // we have to rotate the whole system, so that the old x-axis points towards the new one

    vector<double>& x = c1;
    vector<double>& y = c2;
    vector<double>& z = c3;
    

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

CylindricalCoordinateSet::CylindricalCoordinateSet() {
    m_representation = CYLINDRICAL;
    m_size = 0;
}

CylindricalCoordinateSet::CylindricalCoordinateSet(CartesianCoordinateSet& cs) {
    m_representation = CYLINDRICAL;
    m_size = cs.size();
    for(size_t i = 0; i < m_size; ++i)
    {
        double x, y, z;
        CartesianCoor3D cac(cs.c1[i],cs.c2[i],cs.c3[i]);
        CylinderCoor3D cyc(cac);
                
        c1.push_back(cyc.r);
        c2.push_back(cyc.phi);
        c3.push_back(cyc.z);
    }
}

SphericalCoordinateSet::SphericalCoordinateSet() {
    m_representation = SPHERICAL;
    m_size = 0;    
}

SphericalCoordinateSet::SphericalCoordinateSet(CartesianCoordinateSet& cs) {
    m_representation = SPHERICAL;
    m_size = cs.size();
    for(size_t i = 0; i < m_size; ++i)
    {
        double x, y, z;
        CartesianCoor3D cac(cs.c1[i],cs.c2[i],cs.c3[i]);
        SphericalCoor3D csc(cac);
                
        c1.push_back(csc.r);
        c2.push_back(csc.phi);
        c3.push_back(csc.theta);
    }
}



// end of file
