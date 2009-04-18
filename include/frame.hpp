/*
 *  frame.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef FRAME_HPP_
#define FRAME_HPP_

// common header
#include "common.hpp"

// standard header
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// other headers
#include "atoms.hpp"
#include "coor3d.hpp"

class DcdHeader {	
public:
	int32_t headsize;
	int32_t fingerprint;
	int32_t number_of_frames;
	int32_t dummy1;
	int32_t timesteps_between_frames;
	char buf1[24];
	float size_of_timestep;
	int32_t flag_ext_block1;
	int32_t flag_ext_block2;
	// size1 == size2
};

//forward declaration...
class Atoms;
class Atomselection;


class DcdFrame {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & block1;
        ar & block2;
        ar & x;
        ar & y;
        ar & z;
        ar & number_of_atoms;
        ar & unit_cell_status;
        ar & origin;
    }
	///////////////////
public:
	// block1 contains a triangular matrix:
	// X1 Y1 Y2 Z1 Z2 Z3
	std::vector<double> block1;
	std::vector<double> block2;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	boost::numeric::ublas::matrix<double> coord3Dmatrix;

	int number_of_atoms;
	
	bool unit_cell_status;

	DcdFrame() {}
	~DcdFrame() {}

	std::vector<CartesianCoor3D> unit_cell();
	CartesianCoor3D origin;
	
	double unit_cell_volume() { std::vector<CartesianCoor3D> uc = unit_cell(); return (uc[0].cross_product(uc[1]))*uc[2]; }
	bool has_unit_cell() { return unit_cell_status; }

	void clear() { block1.clear(); block2.clear(); x.clear(); y.clear(); z.clear(); coord3Dmatrix.clear(); unit_cell_status=false;}

//	cartcoord_t coord(int i) { return cartcoord_t(x[i],y[i],z[i]); }
	CartesianCoor3D coord3D(int i); 
	CartesianCoor3D cofm(Atoms& atoms, Atomselection& as); 
	void wrap();
//	void shift();		
};


#endif
