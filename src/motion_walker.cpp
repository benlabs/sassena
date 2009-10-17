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
#include "motion_walker.hpp"
#include "vector_unfold.hpp"
#include "log.hpp"

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/variate_generator.hpp>

using namespace std;

// Brownian Motion

BrownianMotionWalker::BrownianMotionWalker(double displace,long seed,CartesianCoor3D direction) {
	m_seed = seed;
	m_direction = direction;
	m_displace = displace;
	
	p_svu = new SphereVectorUnfold(direction/direction.length());
	p_svu->set_resolution(1000000);
	p_svu->set_seed(m_seed);
	p_svu->set_vectors("mcboostunisphere");
	p_svu->execute();
	
	boost::mt19937 brownian_displace_rng; // that's my random number generator
	brownian_displace_rng.seed(1000000);		
	boost::normal_distribution<double> gauss; // that's my distribution
	p_mynormaldistribution = new boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >(brownian_displace_rng,gauss);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > asdfasd(brownian_displace_rng,gauss); 
}

BrownianMotionWalker::~BrownianMotionWalker() {
	delete p_svu;
    delete p_mynormaldistribution;
}



CartesianCoor3D BrownianMotionWalker::translation(size_t timepos) {
	if (translations.find(timepos)==translations.end()) {
		generate(timepos);
	}
	return translations[timepos];
}


void BrownianMotionWalker::generate(size_t timepos) {
	size_t oldtimepos = translations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	CartesianCoor3D oldtranslation(0,0,0);
	if (oldtimepos>0) {
		oldtranslation = translations[oldtimepos-1];		
	}

	if (timepos>p_svu->vectors().size()) {
		p_svu->set_resolution(p_svu->vectors().size()+timepos+1000000); // this will ensure that we have enough vectors
	}
    
	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
        double normran = (*p_mynormaldistribution)();
        double displacement = m_displace * normran;
		translations[ti]= oldtranslation + displacement*p_svu->vectors()[ti] ;
		oldtranslation = translations[ti];
	}
}


// random walk 

RandomMotionWalker::RandomMotionWalker(double displace,long seed,CartesianCoor3D direction) {
	m_seed = seed;
	m_direction = direction;
	m_displace = displace;
	
	p_svu = new SphereVectorUnfold(direction/direction.length());
	p_svu->set_resolution(1000000);
	p_svu->set_seed(m_seed);
	p_svu->set_vectors("mcboostunisphere");
	p_svu->execute();
}

RandomMotionWalker::~RandomMotionWalker() {
	delete p_svu;
}



CartesianCoor3D RandomMotionWalker::translation(size_t timepos) {
	if (translations.find(timepos)==translations.end()) {
		generate(timepos);
	}
	return translations[timepos];
}


void RandomMotionWalker::generate(size_t timepos) {
	size_t oldtimepos = translations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	CartesianCoor3D oldtranslation(0,0,0);
	if (oldtimepos>0) {
		oldtranslation = translations[oldtimepos-1];		
	}

	if (timepos>p_svu->vectors().size()) {
		p_svu->set_resolution(p_svu->vectors().size()+timepos+1000000); // this will ensure that we have enough vectors
	}

	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
		translations[ti]= oldtranslation + m_displace*p_svu->vectors()[ti] ;
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


