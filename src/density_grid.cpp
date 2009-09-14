/*
 *  density_grid.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "density_grid.hpp"

// standard header
#include <cmath>

// special library headers

// other headers


using namespace std;
	
void DensityGrid::set(CoordinateSet& coordset) {
	for(size_t i = 0; i < coordset.size(); ++i)
	{
		Gridkey3D key = get_cell(CartesianCoor3D(coordset.x[i],coordset.y[i],coordset.z[i]));
		(*this)[key] = true;
	}
}

void DensityGrid::unset(CoordinateSet& coordset) {
	for(size_t i = 0; i < coordset.size(); ++i)
	{
		Gridkey3D key = get_cell(CartesianCoor3D(coordset.x[i],coordset.y[i],coordset.z[i]));
		(*this)[key] = false;
	}
}

void DensityGrid::unset(CoordinateSet& coordset, double nullrange) {
	double a_d = box[0].length()/a_cells;
	double b_d = box[1].length()/b_cells;
	double c_d = box[2].length()/c_cells;

	long a_max = static_cast<long>(nullrange/a_d);
	long a_min = -a_max;
	long b_max = static_cast<long>(nullrange/b_d);
	long b_min = -b_max;
	long c_max = static_cast<long>(nullrange/c_d);
	long c_min = -c_max;
	
	for(size_t ai = a_min; ai <= a_max; ++ai)
	{
		for(size_t bi = b_min; bi <= b_max; ++bi)
		{
			for(size_t ci = c_min; ci <= c_max; ++ci)
			{
				for(size_t i = 0; i < coordset.size(); ++i)
				{		
					CartesianCoor3D c(coordset.x[i]+ai*a_d,coordset.y[i]+bi*b_d,coordset.z[i]+ci*c_d);

					Gridkey3D key = get_cell(c);
					(*this)[key] = false;
				}

			}
		}
	}

}
	
double DensityGrid::coverage() {
	size_t count = 0;
	for (int i=0;i<a_cells;i++) {
		for (int j=0;j<b_cells;j++) {				
			for (int k=0;k<c_cells;k++) {
				if ((*this)[Gridkey3D(i,j,k)]) count++;
			}
		}
	}
	return (count*1.0/(a_cells*b_cells*c_cells));
}

// end of file