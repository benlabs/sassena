/** \file
This file contains a class which defines the contents of a generic frame.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "sample/frame.hpp"

// standard header

// special library headers

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"

using namespace std;

void Frame::clear() { 
	x.clear(); 
	y.clear(); 
	z.clear(); 
	unitcell.clear();
}

CartesianCoor3D Frame::coord3D(size_t i) { 
	return CartesianCoor3D(x[i],y[i],z[i]);	
}

// end of file

