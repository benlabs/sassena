/*
 *  density_grid.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef DENSITY_GRID_HPP_
#define DENSITY_GRID_HPP_

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
#include "atomselection.hpp"
#include "coordinate_set.hpp"
#include "coor3d.hpp"
#include "grid.hpp"


class DensityGrid : public Grid3D<bool> {
public:
	DensityGrid(std::vector<CartesianCoor3D>& uc,double spacing, CartesianCoor3D o, bool def) : Grid3D<bool>(uc,spacing,o,def) {};
	void set(CoordinateSet& coordset);
	void unset(CoordinateSet& coordset);
	void unset(CoordinateSet& coordset,double nullrange);
	
	
	double coverage();
};

#endif 

// end of file