/*
 *  vector_unfold.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef VECTOR_UNFOLD_HPP_
#define VECTOR_UNFOLD_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "coor3d.hpp"

class VectorUnfold {
protected:
	CartesianCoor3D m_q;
	std::vector<CartesianCoor3D> m_unfolded_vectors;	
public:
	virtual void execute() = 0;
	virtual std::vector<CartesianCoor3D>& vectors() = 0;
	
	void clear() { m_unfolded_vectors.clear(); }	
};

class NoVectorUnfold : public VectorUnfold {
public:
	NoVectorUnfold(CartesianCoor3D q);
	
	void execute();
	std::vector<CartesianCoor3D>& vectors();	
};

class SphereVectorUnfold : public VectorUnfold {
	uint32_t m_seed;
	size_t m_resolution;
	std::string m_vectors;
public:
	SphereVectorUnfold(CartesianCoor3D q);
	
	void execute();
	std::vector<CartesianCoor3D>& vectors();	
	
	void set_resolution(size_t resolution);
	void set_seed(uint32_t seed);
	void set_vectors(std::string vectors);
};

class CylinderVectorUnfold : public VectorUnfold {
	CartesianCoor3D m_axis;
	uint32_t m_seed;
	size_t m_resolution;
	std::string m_vectors;
public:
	CylinderVectorUnfold(CartesianCoor3D q);
	
	void execute();
	std::vector<CartesianCoor3D>& vectors();
	
	void set_resolution(size_t resolution);
	void set_seed(uint32_t seed);
	void set_vectors(std::string vectors);
	void set_axis(CartesianCoor3D axis);
};

// we can actually support tons of different unfolding schemes:
class GaussianVectorUnfold : public VectorUnfold {
	
};

// this reads vectors from a file. Only their direction is used to generate the final vectors
class FileVectorUnfold : public VectorUnfold {
public:
	FileVectorUnfold(CartesianCoor3D q);
	
	void execute();
	std::vector<CartesianCoor3D>& vectors();
	
};


#endif

// end of file
