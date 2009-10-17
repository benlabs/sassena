/*
 *  motion_walker.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef MOTION_WALKER_HPP_
#define MOTION_WALKER_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/variate_generator.hpp>

// other headers
#include "coor3d.hpp"
#include "vector_unfold.hpp"

class MotionWalker {

public:
	
	virtual CartesianCoor3D translation(size_t timepos) = 0;	
};

//

class LinearMotionWalker : public MotionWalker {
	CartesianCoor3D m_translate;
	
	void generate(size_t timepos);
public:
	LinearMotionWalker(double displace,CartesianCoor3D direction);
	
	CartesianCoor3D translation(size_t timepos);
};

class FixedMotionWalker : public MotionWalker {
	CartesianCoor3D m_translate;
public:
	FixedMotionWalker(double displace,CartesianCoor3D direction);
	
	CartesianCoor3D translation(size_t timepos);
};

class OscillationMotionWalker : public MotionWalker {

	CartesianCoor3D m_translate;
	double m_frequency;

public:
	OscillationMotionWalker(double displace,double frequency, CartesianCoor3D direction);
	
	CartesianCoor3D translation(size_t timepos);
};


class RandomMotionWalker : public MotionWalker {

	std::map<size_t,CartesianCoor3D> translations;
	
	SphereVectorUnfold* p_svu;
	
	double m_displace;
	long m_seed;
	
	CartesianCoor3D m_direction;
	
	void generate(size_t timepos);
public:
	RandomMotionWalker(double displace,long seed, CartesianCoor3D direction);
	~RandomMotionWalker();
	
	
	CartesianCoor3D translation(size_t timepos);
};


class BrownianMotionWalker : public MotionWalker {

	std::map<size_t,CartesianCoor3D> translations;
	
	SphereVectorUnfold* p_svu;
	
	double m_displace;
	long m_seed;
	
	CartesianCoor3D m_direction;
    boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >* p_mynormaldistribution;
	
	void generate(size_t timepos);
public:
	BrownianMotionWalker(double displace,long seed, CartesianCoor3D direction);
	~BrownianMotionWalker();
	
	
	CartesianCoor3D translation(size_t timepos);
};

#endif

// end of file
