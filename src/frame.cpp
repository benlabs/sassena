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

void Frame::clear() { 
	x.clear(); 
	y.clear(); 
	z.clear(); 
	unitcell.clear();
}


void Frame::push_selection(Atomselection& as) {
	// create empty coordinate set
	CoordinateSet& cs = coordinate_sets[as.name];

	// try to take advantage of hardware pipelining, each array individually
	for (size_t i=0;i<number_of_atoms;i++) {
			if (as.booleanarray[i]) cs.x.push_back(x[i]);
	}

	for (size_t i=0;i<number_of_atoms;i++) {
			if (as.booleanarray[i]) cs.y.push_back(y[i]);
	}

	for (size_t i=0;i<number_of_atoms;i++) {
			if (as.booleanarray[i]) cs.z.push_back(z[i]);
	}

	for (size_t i=0;i<number_of_atoms;i++) {
			if (as.booleanarray[i]) cs.indexes.push_back(i);
	}
	
	
}

void Frame::push_selections(std::vector<Atomselection>& ass) {

	// temp asssociation table
	std::vector<pair<CoordinateSet*,Atomselection*> > cooras;

	for (std::vector<Atomselection>::iterator asi=ass.begin();asi!=ass.end();asi++) {
		// create empty coordinate set
		CoordinateSet* csp = &(coordinate_sets[asi->name]);
		cooras.push_back(pair<CoordinateSet*,Atomselection*>(csp,&(*asi)));
	}
	
	// try to take advantage of hardware pipelining, each array individually
	for (size_t i=0;i<number_of_atoms;i++) {
		for (std::vector<pair<CoordinateSet*,Atomselection*> >::iterator csi=cooras.begin();csi!=cooras.end();csi++) {
			if (csi->second->booleanarray[i]) csi->first->x.push_back(x[i]);
		}
	}

	for (size_t i=0;i<number_of_atoms;i++) {
		for (std::vector<pair<CoordinateSet*,Atomselection*> >::iterator csi=cooras.begin();csi!=cooras.end();csi++) {
			if (csi->second->booleanarray[i]) csi->first->y.push_back(y[i]);
		}
	}

	for (size_t i=0;i<number_of_atoms;i++) {
		for (std::vector<pair<CoordinateSet*,Atomselection*> >::iterator csi=cooras.begin();csi!=cooras.end();csi++) {
			if (csi->second->booleanarray[i]) csi->first->z.push_back(z[i]);
		}
	}

	for (size_t i=0;i<number_of_atoms;i++) {
		for (std::vector<pair<CoordinateSet*,Atomselection*> >::iterator csi=cooras.begin();csi!=cooras.end();csi++) {
			if (csi->second->booleanarray[i]) csi->first->indexes.push_back(i);
		}
	}
	
}

CartesianCoor3D Frame::coord3D(size_t i) { 
	return CartesianCoor3D(x[i],y[i],z[i]);	
}

CartesianCoor3D Frame::cofm(Atoms& atoms, Atomselection& as) {

	if (as.empty()) {
		cerr << "Warning! Computing Center of Mass for an empty atomselection" << endl;
		cerr << "Setting Center of mass to (0,0,0)" << endl;		
		return CartesianCoor3D(0,0,0);
	}
	

	double xt,yt,zt,m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	for (size_t i=0;i<as.size();i++) {
		mi = atoms[as[i]].mass;
		m += mi;
		xt += x[as[i]]*mi;
		yt += y[as[i]]*mi;
		zt += z[as[i]]*mi;
	}
	
	return CartesianCoor3D(xt/m,yt/m,zt/m);
}

// Map each atom back to the origin unit cell
void Frame::wrap() {
	std::vector<CartesianCoor3D>& uc = unitcell;
	double uc0l = uc[0].length();
	double uc1l = uc[1].length();
	double uc2l = uc[2].length();
	CartesianCoor3D uc0n = uc[0]/uc0l;
	CartesianCoor3D uc1n = uc[1]/uc1l;
	CartesianCoor3D uc2n = uc[2]/uc2l;

	for (size_t i=0;i<number_of_atoms;i++) {
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


// end of file

