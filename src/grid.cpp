/*
 *  grid.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "grid.hpp"

// standard header
#include <cmath>

// special library headers

// other headers


using namespace std;

bool Gridkey3D::operator<(const Gridkey3D that) const {
	if (this->ix==that.ix) {
		if (this->iy==that.iy) {
			if (this->iz==that.iz) return false; else return this->iz < that.iz;
		}
		return this->iy < that.iy;
	}
	return this->ix < that.ix;
	
}


//template<class V> Grid3D<V>::Grid3D(int r,vector<CartesianCoor3D>& uc,CartesianCoor3D o) {
//	res = r;
//	box = uc;
//	origin = o;
//	for (int i=0;i<res;i++) {
//		for (int j=0;j<res;j++) {				
//			for (int k=0;k<res;k++) {
//				this->operator[](Gridkey3D(i,j,k)); // initialize with empty value type
//			}
//		}
//	}
//}
//
//
//template<class V> Grid3D<V>::Grid3D(int r,vector<CartesianCoor3D>& uc,CartesianCoor3D o, V v) {
//	res = r;
//	box = uc;
//	origin = o;
//	for (int i=0;i<res;i++) {
//		for (int j=0;j<res;j++) {				
//			for (int k=0;k<res;k++) {
//				this->operator[](Gridkey3D(i,j,k)) = v; // initialize with value type
//			}
//		}
//	}
//}


template<class V> Grid3D<V>::Grid3D(vector<CartesianCoor3D>& uc,double spacing,CartesianCoor3D o, V v) {
	box = uc;
	origin = o;
	slice_box(spacing); // now a_cells, b_cells, c_cells is set
	for (int i=0;i<a_cells;i++) {
		for (int j=0;j<b_cells;j++) {				
			for (int k=0;k<c_cells;k++) {
				this->operator[](Gridkey3D(i,j,k)) = v; // initialize with value type
			}
		}
	}
}

template<class V> void Grid3D<V>::slice_box(double spacing) {
	double a_cells = ((box[0]).length() / spacing); // then box[0].length / a_cells is the real spacing!
	double b_cells = ((box[1]).length() / spacing); // then box[0].length / a_cells is the real spacing!
	double c_cells = ((box[2]).length() / spacing); // then box[0].length / a_cells is the real spacing!	
}

// determine gridkey element given a cartesian Coordinate...
// it is illegal to call this function with a c being not inside the defined box (box,origin)
template<class V> Gridkey3D Grid3D<V>::get_cell(CartesianCoor3D c) {
	CartesianCoor3D uc0 = box[0];
	CartesianCoor3D uc1 = box[1];
	CartesianCoor3D uc2 = box[2];
	double uc0l = uc0.length();
	double uc1l = uc1.length();
	double uc2l = uc2.length();
	CartesianCoor3D uc0n = uc0/uc0l;
	CartesianCoor3D uc1n = uc1/uc1l;
	CartesianCoor3D uc2n = uc2/uc2l;
	CartesianCoor3D uc0d = uc0/a_cells;
	CartesianCoor3D uc1d = uc1/b_cells;
	CartesianCoor3D uc2d = uc2/c_cells;
	
	CartesianCoor3D neworigin = -0.5*uc0 - 0.5*uc1 - 0.5*uc2;
	CartesianCoor3D cold = c;
 	c = c - neworigin;
//	c = c - 0.5*uc0;
//	c = c - 0.5*uc1;
//	c = c - 0.5*uc2;
//	cout <<  "atom now has: " << c << " ";
			
	double cuc0n = (c*uc0n);
	double cuc1n = (c*uc1n);
	double cuc2n = (c*uc2n);		


	int iu0 = int(cuc0n/uc0d.length());
	int iu1 = int(cuc1n/uc1d.length());
	int iu2 = int(cuc2n/uc2d.length());
	
	// rare case: atom sits on the "right" border...
//	if ((iu0==resolution()) && (cuc0n==(uc0d.length()*resolution()))) iu0--;
//	if ((iu1==resolution()) && (cuc1n==(uc1d.length()*resolution()))) iu1--;
//	if ((iu2==resolution()) && (cuc2n==(uc2d.length()*resolution()))) iu2--;		
		
	
//	cout << "iu0=" << iu0 << endl;
//	cout << "iu1=" << iu1 << endl;
//	cout << "iu2=" << iu2 << endl;
	
	if ( (iu0<0) || (iu1<0) || (iu2<0) || (iu0>=a_cells) || (iu1>=b_cells) || (iu2>=c_cells) ) {
		cerr << "grid index out of range" << endl;
		cerr << "Atom coordinates: " << cold << endl;		
		cerr << "Atom coordinates: " << c << endl;
		cerr << "indices: " << iu0 << ", " << iu1 << ", " << iu2 << endl;
		cerr << "origin:" << origin << endl;
		cerr <<  cuc0n/uc0d.length() << endl;
		cerr <<  cuc1n/uc1d.length() << endl;
		cerr <<  cuc2n/uc2d.length() << endl;
		cerr <<  "a_cells: "<< a_cells << endl;
		cerr <<  "b_cells: "<< b_cells << endl;
		cerr <<  "c_cells: "<< c_cells << endl;
		throw("");		
	}
	
	return Gridkey3D(iu0,iu1,iu2);
}	


// determine gridkey elements given a origin and a box with ...
template<class V> vector<Gridkey3D> Grid3D<V>::get_cells(CartesianCoor3D c,CartesianCoor3D sc) {
	vector<Gridkey3D> gks;
	
	// assume linear relationship between grid coordinates and cartesian coordinates....
	
	CartesianCoor3D sc0 = CartesianCoor3D(sc.x,0,0);
	CartesianCoor3D sc1 = CartesianCoor3D(0,sc.y,0);
	CartesianCoor3D sc2 = CartesianCoor3D(0,0,sc.z);
	
	CartesianCoor3D uc0 = box[0];
	CartesianCoor3D uc1 = box[1];
	CartesianCoor3D uc2 = box[2];
	double uc0l = uc0.length();
	double uc1l = uc1.length();
	double uc2l = uc2.length();
	CartesianCoor3D uc0n = uc0/uc0l;
	CartesianCoor3D uc1n = uc1/uc1l;
	CartesianCoor3D uc2n = uc2/uc2l;
	CartesianCoor3D uc0d = uc0/a_cells;
	CartesianCoor3D uc1d = uc1/b_cells;
	CartesianCoor3D uc2d = uc2/c_cells;


	//	cout <<  "atom with " << c << "has ";

	CartesianCoor3D neworigin = -0.5*uc0 - 0.5*uc1 - 0.5*uc2;
 	c = c - neworigin;
	//	c = c - 0.5*uc0;
	//	c = c - 0.5*uc1;
	//	c = c - 0.5*uc2;
	//	cout <<  "atom now has: " << c << " ";

	CartesianCoor3D c0_l = c - sc0;
	CartesianCoor3D c0_r = c + sc0;
	CartesianCoor3D c1_l = c - sc1;
	CartesianCoor3D c1_r = c + sc1;
	CartesianCoor3D c2_l = c - sc2;
	CartesianCoor3D c2_r = c + sc2;

	int iu0_l = int((c0_l*uc0n)/uc0d.length());
	int iu0_r = int((c0_r*uc0n)/uc0d.length());
	int iu1_l = int((c1_l*uc1n)/uc1d.length());
	int iu1_r = int((c1_r*uc1n)/uc1d.length());
	int iu2_l = int((c2_l*uc2n)/uc2d.length());
	int iu2_r = int((c2_r*uc2n)/uc2d.length());
	
	for (int i0=iu0_l;i0<=iu0_r;i0++) {
		for (int i1=iu1_l;i1<=iu1_r;i1++) {
			for (int i2=iu2_l;i2<=iu2_r;i2++) {
				if ( (i0>=0) && (i0<a_cells) && 
					 (i1>=0) && (i1<b_cells) && 
					 (i2>=0) && (i2<c_cells) )
					gks.push_back(Gridkey3D(i0,i1,i2));
			}
		}		
	}
	
	return gks;
}	


template<class V> std::ostream& operator<<(std::ostream& os, Grid3D<V>& grid) {
	for (int i=0;i<grid.a_cells;i++) {
		for (int j=0;j<grid.b_cells;j++) {				
			for (int k=0;k<grid.c_cells;k++) {
				os << grid[Gridkey3D(i,j,k)] << endl;
			}
		}
	}
	return os;
}


 //special output for vectorized grids..
ostream& operator<<(ostream& os, vGrid3D& grid) {
	int totalcount =0;
	for (int i=0;i<grid.a_cells;i++) {
		for (int j=0;j<grid.b_cells;j++) {				
			for (int k=0;k<grid.c_cells;k++) {
				os << "(" << i << "," << j << "," << k << ")= ";
				vector<int>& v = grid[Gridkey3D(i,j,k)];
				int count=0;
				for (vector<int>::iterator vi=v.begin();vi!=v.end();vi++) {
//					os << *vi;
					count++; totalcount++;
				}
				os << count << " atoms";
				os << endl;
			}
		}
	}
	os << "total: " << totalcount << " atoms" << endl;
	return os;
}

// explicit instantiation
template class Grid3D<vector<int> >;
template class Grid3D<bool >;

// end of file
