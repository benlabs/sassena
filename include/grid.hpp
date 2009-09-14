/*
 *  grid.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef GRID_HPP_
#define GRID_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>

// special library headers
#include <libconfig.h++>
#include <boost/regex.hpp>

// other headers
#include "coor3d.hpp"

class Gridkey3D {
public:
	int ix,iy,iz;
	
	Gridkey3D(int a,int b,int c) : ix(a), iy(b), iz(c) {} 
	
	bool operator<(const Gridkey3D that) const;
};

template<class V> class Grid3D : public std::map< Gridkey3D, V > {
protected:
	std::vector<CartesianCoor3D> box;
	CartesianCoor3D origin;
	
	void slice_box(double spacing);
	
public:
		size_t a_cells,b_cells,c_cells;
	Grid3D(std::vector<CartesianCoor3D>& uc,double spacing, CartesianCoor3D o, V v);
//	Grid3D(std::vector<CartesianCoor3D>& uc,double spacing, CartesianCoor3D o);
	
//	int resolution();
	double total_volume() { return (box[0].cross_product(box[1]))*box[2]; }
	double element_volume() { return total_volume()/(a_cells*b_cells*c_cells); }
	
	Gridkey3D get_cell(CartesianCoor3D c);
	std::vector<Gridkey3D> get_cells(CartesianCoor3D c,CartesianCoor3D sc);		
	//specialization:
	//std::ostream& operator<<(std::ostream& os);
};

template<class V> std::ostream& operator<<(std::ostream& os, Grid3D<V>& grid);

// handy shortcuts:
typedef Grid3D<std::vector<int> > vGrid3D;
std::ostream& operator<<(std::ostream& os, vGrid3D& grid);

#endif
