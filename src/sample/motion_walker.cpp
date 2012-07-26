/** \file
This file contains a set of artifical motion generators. They are used to superimpose artifical with the specific dynamics provided by the trajectory files.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "sample/motion_walker.hpp"
#include "control.hpp"
#include "log.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	

#include <boost/random/variate_generator.hpp>

using namespace std;



// Angular Brownian Motion

RotationalBrownianMotionWalker::RotationalBrownianMotionWalker(double displace,unsigned long seed,long sampling): m_init(true) {
	m_seed = seed;
    m_sampling = sampling;
	m_displace = displace;
}

void RotationalBrownianMotionWalker::init() {
	boost::mt19937 brownian_displace_rng; // that's my random number generator
	brownian_displace_rng.seed(m_seed);

	boost::normal_distribution<double> gauss; // that's my distribution

	p_mynormaldistribution = new boost::variate_generator<boost::mt19937, boost::normal_distribution<double> >(brownian_displace_rng,gauss);
	
    m_init = false;
}

RotationalBrownianMotionWalker::~RotationalBrownianMotionWalker() {
    if (!m_init) {
        delete p_mynormaldistribution;        
    }
}

boost::numeric::ublas::matrix<double> RotationalBrownianMotionWalker::transform(size_t timepos) {
	if (transformations.find(timepos)==transformations.end()) {
		generate(timepos);
	}
	return transformations[timepos];
}

void RotationalBrownianMotionWalker::generate(size_t timepos) {
    if (m_init) init();    
	size_t oldtimepos = transformations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	boost::numeric::ublas::matrix<double> oldtransform;
	oldtransform=boost::numeric::ublas::identity_matrix<double>(4,4);
	if (oldtimepos>0) {
		oldtransform = transformations[oldtimepos-1];		
	}
    
	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
        double normran1 = (*p_mynormaldistribution)();
        double normran2 = (*p_mynormaldistribution)();
        double normran3 = (*p_mynormaldistribution)();
	    for(size_t i = 0; i < (m_sampling-1); ++i)
        {
            (*p_mynormaldistribution)();
        }
	    
        double angle1 = m_displace * normran1 * M_PI / 180; 
        double angle2 = m_displace * normran2 * M_PI / 180; 
        double angle3 = m_displace * normran3 * M_PI / 180; 

		boost::numeric::ublas::matrix<double> newtransform,r1,r2,r3,r12; 
		r1=r2=r3=newtransform=boost::numeric::ublas::identity_matrix<double>(4,4);
		r1(0,0)=cos(angle1); 	r1(0,1)=sin(angle1);
		r1(1,0)=-sin(angle1);	r1(1,1)=cos(angle1);
		r2(0,0)=cos(angle2); 	r2(0,2)=-sin(angle2);
		r2(2,0)=sin(angle2);	r2(2,2)=cos(angle2);
		r3(1,1)=cos(angle3); 	r3(1,2)=sin(angle3);
		r3(2,1)=-sin(angle3);	r3(2,2)=cos(angle3);
		r12=boost::numeric::ublas::prod(r1,r2);
		newtransform=boost::numeric::ublas::prod(r12,r3);
		transformations[ti]= boost::numeric::ublas::prod(oldtransform,newtransform);
		oldtransform = transformations[ti];
	}
}

// Brownian Motion

BrownianMotionWalker::BrownianMotionWalker(double displace,unsigned long seed,long sampling, CartesianCoor3D direction): m_init(true) {
	m_seed = seed;
    m_sampling = sampling;
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

boost::numeric::ublas::matrix<double> BrownianMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
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
	    for(size_t i = 0; i < (m_sampling-1); ++i)
        {
            (*p_mynormaldistribution)();
            (*p_myspheredistribution)();
        }
	    
        double displacement = m_displace * normran; 
		translations[ti]= oldtranslation + displacement*CartesianCoor3D(sphereran[0],sphereran[1],sphereran[2]) ;
		oldtranslation = translations[ti];
	}
}

// Localized Brownian Motion

LocalBrownianMotionWalker::LocalBrownianMotionWalker(double radius, double displace,unsigned long seed,long sampling, CartesianCoor3D direction): m_init(true) {
	m_seed = seed;
    m_radius = radius;
    m_sampling = sampling;
	m_direction = direction;
	m_displace = displace;
	
	if (m_displace>m_radius) {
        Err::Inst()->write("radius size for local brownian motion smaller than displacement!");
        throw;
    }
}

void LocalBrownianMotionWalker::init() {
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

LocalBrownianMotionWalker::~LocalBrownianMotionWalker() {
    if (!m_init) {
    	delete p_myspheredistribution;
        delete p_mynormaldistribution;        
    }
}

CartesianCoor3D LocalBrownianMotionWalker::translation(size_t timepos) {
    
	if (translations.find(timepos)==translations.end()) {
		generate(timepos);
	}
	return translations[timepos];
}

boost::numeric::ublas::matrix<double> LocalBrownianMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
}

void LocalBrownianMotionWalker::generate(size_t timepos) {
    if (m_init) init();    
	size_t oldtimepos = translations.size();
	if (timepos<oldtimepos) return; // don't generate anything if we already have it

	CartesianCoor3D oldtranslation(0,0,0);
	if (oldtimepos>0) {
		oldtranslation = translations[oldtimepos-1];		
	}
    
	for(size_t ti = oldtimepos; ti <= timepos; ++ti)
	{
        CartesianCoor3D newtrans;
        while(true) {
            double normran = (*p_mynormaldistribution)();
            vector<double> sphereran = (*p_myspheredistribution)();
    	    for(size_t i = 0; i < (m_sampling-1); ++i)
            {
                (*p_mynormaldistribution)();
                (*p_myspheredistribution)();
            }
            double displacement = m_displace * normran; 
    		newtrans = oldtranslation + displacement*CartesianCoor3D(sphereran[0],sphereran[1],sphereran[2]) ;
            if (newtrans.length()>m_radius) continue; else break;
        }
        translations[ti] = newtrans;
        oldtranslation = newtrans;
	}
}

// Localized Brownian Motion


// random walk 

RandomMotionWalker::RandomMotionWalker(double displace,unsigned long seed,long sampling, CartesianCoor3D direction): m_init(true) {
	m_seed = seed;
    m_sampling = sampling;
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

boost::numeric::ublas::matrix<double> RandomMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
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
        for(size_t i = 0; i < (m_sampling-1); ++i)
        {
            vector<double> discard_sphereran = (*p_myspheredistribution)();
        }
        
		translations[ti]= oldtranslation + m_displace*CartesianCoor3D(sphereran[0],sphereran[1],sphereran[2]) ;
		oldtranslation = translations[ti];
	}
}


// oscillation

OscillationMotionWalker::OscillationMotionWalker(double displace,double frequency,long sampling, CartesianCoor3D direction) {
	m_translate = displace*direction/direction.length();
	m_frequency = frequency;
    m_sampling = sampling;
}

CartesianCoor3D OscillationMotionWalker::translation(size_t timepos) {
	return m_translate*sin(2*M_PI*timepos*m_frequency*m_sampling);
}

boost::numeric::ublas::matrix<double> OscillationMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
}

// linear motion

LinearMotionWalker::LinearMotionWalker(double displace,long sampling, CartesianCoor3D direction) {
	m_translate = displace*sampling*direction/direction.length();
}

CartesianCoor3D LinearMotionWalker::translation(size_t timepos) {
	return timepos*m_translate;
}

boost::numeric::ublas::matrix<double> LinearMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
}

// fixed point translation

FixedMotionWalker::FixedMotionWalker(double displace,CartesianCoor3D direction) {
	m_translate = displace*direction/direction.length();
}

CartesianCoor3D FixedMotionWalker::translation(size_t timepos) {
	return m_translate;
}

boost::numeric::ublas::matrix<double> FixedMotionWalker::transform(size_t timepos) {

    using namespace boost::numeric::ublas;
	matrix<double> T(4,4); T= identity_matrix<double>(4);
	CartesianCoor3D translation = this->translation(timepos);
	T(3,0)=translation.x;
	T(3,1)=translation.y;
	T(3,2)=translation.z;
	return T;
}

// end of file
