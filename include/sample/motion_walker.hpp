/** \file
This file contains a set of artifical motion generators. They are used to superimpose artifical with the specific dynamics provided by the trajectory files.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>

// other headers
#include "math/coor3d.hpp"

/** 
Interface for an artficial motion generator
*/
class MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        
    }

public:
	
	virtual boost::numeric::ublas::matrix<double> transform(size_t timepos) = 0;	
};

/** 
Generates a linear motion
*/
class LinearMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<MotionWalker, LinearMotionWalker>(*this);
        ar & m_translate;
    }    
	CartesianCoor3D m_translate;
	
	void generate(size_t timepos);
    LinearMotionWalker() {}
	CartesianCoor3D translation(size_t timepos);

public:
	LinearMotionWalker(double displace,long sampling, CartesianCoor3D direction);

	boost::numeric::ublas::matrix<double> transform(size_t timepos);	
};

/** 
No motion. Fixed offset.
*/
class FixedMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<MotionWalker, FixedMotionWalker>(*this);
        ar & m_translate;        
    }    
    
	CartesianCoor3D m_translate;
    FixedMotionWalker() {}	

	CartesianCoor3D translation(size_t timepos);    
public:
	FixedMotionWalker(double displace,CartesianCoor3D direction);
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);	
};

/** 
Sinoidal Oscillation with amplitude and frequency.
*/
class OscillationMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<MotionWalker, OscillationMotionWalker>(*this);

        ar & m_translate;        
        ar & m_frequency;
        ar & m_sampling;
    }
	CartesianCoor3D m_translate;
	double m_frequency;
    size_t m_sampling;

	OscillationMotionWalker() {}
	CartesianCoor3D translation(size_t timepos);

public:
	OscillationMotionWalker(double displace,double frequency,long sampling, CartesianCoor3D direction);
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);	
};

/** 
Simple random motion in 3D with a constant step size
*/
class RandomMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        translations.clear(); // don't transmit cache
        ar & boost::serialization::base_object<MotionWalker, RandomMotionWalker>(*this);
        ar & m_displace;
        ar & m_seed;
        ar & m_direction;
        ar & m_sampling;
    }
	std::map<size_t,CartesianCoor3D> translations;
	
    bool m_init; // init flag 
	
	double m_displace;
	unsigned long m_seed;
    size_t m_sampling;
    	
	CartesianCoor3D m_direction;
    boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >* p_myspheredistribution;
    void init();

	void generate(size_t timepos);
	RandomMotionWalker(): m_init(true) {}    

	CartesianCoor3D translation(size_t timepos);

public:
	RandomMotionWalker(double displace,unsigned long seed,long sampling,  CartesianCoor3D direction);
	~RandomMotionWalker();
	
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);
};

/** 
Brownian motion in 3D which draws the step size from a gaussian distribution for the given mean (displace)
*/
class BrownianMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        translations.clear();
        ar & boost::serialization::base_object<MotionWalker, BrownianMotionWalker>(*this);
        ar & m_displace;
        ar & m_seed;
        ar & m_sampling;
        ar & m_direction;        
    }
	std::map<size_t,CartesianCoor3D> translations;

    bool m_init; // init flag 
    
	double m_displace;
	unsigned long m_seed;
    size_t m_sampling;
	
	CartesianCoor3D m_direction;
    boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >* p_mynormaldistribution;
    boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >* p_myspheredistribution;
	
    void init();
    
	void generate(size_t timepos);

    BrownianMotionWalker() : m_init(true) {}	
	CartesianCoor3D translation(size_t timepos);

public:
	
	BrownianMotionWalker(double displace,unsigned long seed,long sampling,  CartesianCoor3D direction);
	~BrownianMotionWalker();
	
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);
};


/** 
Brownian motion in 3D which draws the step size from a gaussian distribution for the given mean (displace) and constrains the motion to the maximum radius.
*/
class LocalBrownianMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        translations.clear();
        ar & boost::serialization::base_object<MotionWalker, LocalBrownianMotionWalker>(*this);
        ar & m_displace;
        ar & m_seed;
        ar & m_radius;
        ar & m_sampling;
        ar & m_direction;        
    }
	std::map<size_t,CartesianCoor3D> translations;

    bool m_init; // init flag 
    
    double m_radius;
	double m_displace;
	unsigned long m_seed;
    size_t m_sampling;
	
	CartesianCoor3D m_direction;
    boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >* p_mynormaldistribution;
    boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >* p_myspheredistribution;
	
    void init();
    
	void generate(size_t timepos);

    LocalBrownianMotionWalker() : m_init(true) {}	
	CartesianCoor3D translation(size_t timepos);

public:
	
	LocalBrownianMotionWalker(double radius, double displace,unsigned long seed,long sampling,  CartesianCoor3D direction);
	~LocalBrownianMotionWalker();
	
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);
};


/** 
Rotational Brownian motion in 3D which draws the angular step size from a gaussian distribution for the given mean (displace)
*/
class RotationalBrownianMotionWalker : public MotionWalker {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        transformations.clear();
        ar & boost::serialization::base_object<MotionWalker, RotationalBrownianMotionWalker>(*this);
        ar & m_displace;
        ar & m_seed;
        ar & m_sampling;
    }
	std::map<size_t,boost::numeric::ublas::matrix<double> > transformations;

    bool m_init; // init flag 
    
	double m_displace;
	unsigned long m_seed;
    size_t m_sampling;
	
    boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >* p_mynormaldistribution;
    boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >* p_myspheredistribution;
	
    void init();
    
	void generate(size_t timepos);

    RotationalBrownianMotionWalker() : m_init(true) {}

public:
	
	RotationalBrownianMotionWalker(double displace,unsigned long seed,long sampling);
	~RotationalBrownianMotionWalker();
	
	
	boost::numeric::ublas::matrix<double> transform(size_t timepos);
};
#endif

// end of file
