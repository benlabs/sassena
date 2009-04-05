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
#include "frame.hpp"

// standard header

// special library headers

// other headers
#include "atoms.hpp"
#include "atomselection.hpp"
#include "coor3d.hpp"


using namespace std;

vector<CartesianCoor3D> DcdFrame::unit_cell() {
	vector<CartesianCoor3D> c;
	c.push_back( CartesianCoor3D(block1[0],0,0)                 );
	c.push_back( CartesianCoor3D(block1[1],block1[2],0)         );
	c.push_back( CartesianCoor3D(block1[3],block1[4],block1[5]) );
	return c;
}

CartesianCoor3D DcdFrame::coord3D(int i) { 
	return CartesianCoor3D(x[i],y[i],z[i]);	
}


CartesianCoor3D DcdFrame::cofm(Atoms& atoms, Atomselection& as) {
	double xt,yt,zt,m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	if (as.empty()) {
		cerr << "Warning! Computing Center of Mass for an empty atomselection" << endl;
		cerr << "Setting Center of mass to (0,0,0)" << endl;		
		return CartesianCoor3D(0,0,0);
	}

	for (Atomselection::iterator asi=as.begin();asi!=as.end();asi++) {
		mi = atoms[*asi].mass;
		m += mi; 
		xt += x[*asi]*mi;	yt += y[*asi]*mi; zt += z[*asi]*mi;
	}


	return CartesianCoor3D(xt/m,yt/m,zt/m);
}

// Map each atom back to the origin unit cell
void DcdFrame::wrap() {
	vector<CartesianCoor3D> uc = unit_cell();
	double uc0l = uc[0].length();
	double uc1l = uc[1].length();
	double uc2l = uc[2].length();
	CartesianCoor3D uc0n = uc[0]/uc0l;
	CartesianCoor3D uc1n = uc[1]/uc1l;
	CartesianCoor3D uc2n = uc[2]/uc2l;

	for (int i=0;i<number_of_atoms;i++) {
		CartesianCoor3D c = coord3D(i) - origin;

	double cuc0n = (c*uc0n);
	double cuc1n = (c*uc1n);
	double cuc2n = (c*uc2n);		
	
	if (cuc0n<0) 
		c = c - int(cuc0n/uc0l-0.5) * uc[0];
	else
		c = c - int(cuc0n/uc0l+0.5) * uc[0];
	if (cuc1n<0)
		c = c - int(cuc1n/uc1l-0.5) * uc[1];
	else
		c = c - int(cuc1n/uc1l+0.5) * uc[1];
	if (cuc2n<0)
		c = c - int(cuc2n/uc2l-0.5) * uc[2];			
	else
		c = c - int(cuc2n/uc2l+0.5) * uc[2];			
		
		if ( (c.x>uc0l/2) || 
		(c.y>uc1l/2) ||
		(c.z>uc2l/2) ||  
		(c.x<-uc0l/2) || 
		(c.y<-uc1l/2) ||
		(c.z<-uc2l/2)) {
			cerr << "wrapping didn't work" << endl;
			cerr << "coordinates:" << x[i] << ", " << y[i] << ", " << z[i] << endl;
			cerr << "coordinates/c:" << c.x << ", " << c.y << ", " << c.z << endl;			
			cerr << "index: " << i<< endl;
			throw("");
		}

	x[i] = c.x;
	y[i] = c.y;
	z[i] = c.z;
		
//	x[i] = c.x + origin.x;
//	y[i] = c.y + origin.y;
//	z[i] = c.z + origin.z;
		
	}
	origin = CartesianCoor3D(0,0,0);
}

// shift system origin to center of mass...
//void DcdFrame::shift() {
//	for (int i=0;i<number_of_atoms;i++) {
//		CartesianCoor3D c = coord3D(i) - origin;
//
//		x[i] = c.x;
//		y[i] = c.y;
//		z[i] = c.z;
//	}
//	origin = CartesianCoor3D(0,0,0);
//}

// end of file

