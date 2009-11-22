/*
 *  motion_walker.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "sample/motion_walker.hpp"
#include "control.hpp"

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	

#include <boost/random/variate_generator.hpp>

using namespace std;

// Brownian Motion

BrownianMotionWalker::BrownianMotionWalker(double displace,long seed,CartesianCoor3D direction): m_init(true) {
	m_seed = seed;
	m_direction = direction;
	m_displace = displace;
}

void BrownianMotionWalker::init() {
	boost::mt19937 brownian_displace_rng; // that's my random number generator
	boost::mt19937 spherical_rng; // that's my random number generator
	brownian_displace_rng.seed(m_seed);
    spherical_rng.seed(m_seed+1);

	boost::normal_distribution<double> gauss; // that's my distribution
	boost::uniform_on_sphere<double> sphere(3); // that's my distribution

	p_mynormaldistribution = new boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >(brownian_displace_rng,gauss);
    p_myspheredistribution = new boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >(spherical_rng,sphere);
	
    m_init = false;
}

BrownianMotionWalker::~BrownianMotionWalker() {
    if (!m_init) {
    	delete p_myspheredistribution;
        delete p_mynormaldistribution;        
    }
}

CartesianCoor3D BrownianMotionWalker::translation(size_t timepos) {
    
	if (translations.find(timepos)==translations.end()) {
		generate(timepos);
	}
	return translations[timepos];
}


void BrownianMotionWalker::generate(size_t timepos) {
    if (m_init) init();    
	size_t oldtimepos = translations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	CartesianCoor3D oldtranslation(0,0,0);
	if (oldtimepos>0) {
		oldtranslation = translations[oldtimepos-1];		
	}
    
	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
        double normran = (*p_mynormaldistribution)();
        vector<double> sphereran = (*p_myspheredistribution)();
        double displacement = m_displace * normran; 
		translations[ti]= oldtranslation + displacement*CartesianCoor3D(sphereran[0],sphereran[1],sphereran[2]) ;
		oldtranslation = translations[ti];
	}
}


// random walk 

RandomMotionWalker::RandomMotionWalker(double displace,long seed,CartesianCoor3D direction): m_init(true) {
	m_seed = seed;
	m_direction = direction;
	m_displace = displace;
}

void RandomMotionWalker::init() {
	
	boost::mt19937 spherical_rng; // that's my random number generator
    spherical_rng.seed(m_seed);

	boost::uniform_on_sphere<double> sphere(3); // that's my distribution
    p_myspheredistribution = new boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> >(spherical_rng,sphere);

    m_init = false;
}


RandomMotionWalker::~RandomMotionWalker() {
    if (!m_init) {    
	    delete p_myspheredistribution;
    }
}



CartesianCoor3D RandomMotionWalker::translation(size_t timepos) {
	if (translations.find(timepos)==translations.end()) {
		generate(timepos);
	}
	return translations[timepos];
}


void RandomMotionWalker::generate(size_t timepos) {
    if (m_init) init();
	size_t oldtimepos = translations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	CartesianCoor3D oldtranslation(0,0,0);
	if (oldtimepos>0) {
		oldtranslation = translations[oldtimepos-1];		
	}

	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
	    vector<double> sphereran = (*p_myspheredistribution)();
		translations[ti]= oldtranslation + m_displace*CartesianCoor3D(sphereran[0],sphereran[1],sphereran[2]) ;
		oldtranslation = translations[ti];
	}
}


// oscillation

OscillationMotionWalker::OscillationMotionWalker(double displace,double frequency,CartesianCoor3D direction) {
	m_translate = displace*direction/direction.length();
	m_frequency = frequency;
}

CartesianCoor3D OscillationMotionWalker::translation(size_t timepos) {
	return m_translate*sin(2*M_PI*timepos*m_frequency);
}


// linear motion

LinearMotionWalker::LinearMotionWalker(double displace,CartesianCoor3D direction) {
	m_translate = displace*direction/direction.length();
}

CartesianCoor3D LinearMotionWalker::translation(size_t timepos) {
	return timepos*m_translate;
}

// fixed point translation

FixedMotionWalker::FixedMotionWalker(double displace,CartesianCoor3D direction) {
	m_translate = displace*direction/direction.length();
}

CartesianCoor3D FixedMotionWalker::translation(size_t timepos) {
	return m_translate;
}


