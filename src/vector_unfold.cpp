/*
 *  vector_unfold.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "vector_unfold.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>


#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>	

// other headers
#include "coor3d.hpp"
#include "geometry.hpp"
#include "log.hpp"
#include "parameters.hpp"

using namespace std;

NoVectorUnfold::NoVectorUnfold(CartesianCoor3D q) {
	m_q = q;
}

void NoVectorUnfold::execute() {
	m_qvectors.push_back(m_q);
}

vector<CartesianCoor3D>& NoVectorUnfold::qvectors() {
	return m_qvectors;
}

FileVectorUnfold::FileVectorUnfold(CartesianCoor3D q) {
	m_q = q;
}

void FileVectorUnfold::execute() {
	m_qvectors.push_back(m_q);
}

vector<CartesianCoor3D>& FileVectorUnfold::qvectors() {
	return m_qvectors;
}

SphereVectorUnfold::SphereVectorUnfold(CartesianCoor3D q) {
	m_q = q;
}

void SphereVectorUnfold::set_resolution(size_t resolution) {
	m_resolution = resolution;
}

void SphereVectorUnfold::set_seed(uint32_t seed) {
	m_seed = seed;
}

void SphereVectorUnfold::set_vectors(string vectors) {
	m_vectors = vectors;
}

vector<CartesianCoor3D>& SphereVectorUnfold::qvectors() {
	return m_qvectors;
}

void SphereVectorUnfold::execute() {
	SphericalCoor3D qs(m_q);
	
	// we can't unfold a point:
	if (qs.r==0.0) {
		m_qvectors.push_back( m_q );
		return;
	}
	
	
	const double M_2PI = 2*M_PI;	
	const double radincr = m_resolution*(M_2PI)/360;		
	const double arcincr = radincr*qs.r;

	int num=0;

	srand(m_seed);

	if (m_vectors=="rasterlinear") {
		for (double t=0;t<M_PI;t+=radincr) {	
			for (double p=0;p<M_2PI;p+=radincr) {
				// take length from scattering vector and angles from unisphere												
				m_qvectors.push_back( SphericalCoor3D(qs.r,p,t) );
		    }
		}		
	}
	else if (m_vectors=="rastersurface") {
		double thetal = M_PI*qs.r; // this is the max distance we "walk on the sphere in theta direction"
		for (double t=0;((t==0) || (t<thetal));t+=arcincr) {	
			double dtheta = (t/thetal)*M_PI;
			double phil = 2*M_PI*sin(dtheta)*qs.r; // this is the max distance we "walk on the sphere in phi direction"		
			if (phil>0) {
				for (double p=0;p<phil;p+=arcincr) {
					double dphi = (p/phil)*2*M_PI;
					// take length from scattering vector and angles from unisphere								
					m_qvectors.push_back( SphericalCoor3D(qs.r,dphi,dtheta) );
			    }
			}
			else {
				m_qvectors.push_back( SphericalCoor3D(qs.r,qs.phi,qs.theta+dtheta) );
			}
		}
	}
	else if (m_vectors=="triangulation") {
		
		DrawSphereHelper d(CartesianCoor3D(0,0,0),qs.r,m_resolution*qs.r);
		d.draw();

		for(vector<CartesianCoor3D>::iterator vi=d.vectors.begin();vi!=d.vectors.end();vi++) {
			// take length from scattering vector and angles from unisphere						
			m_qvectors.push_back( SphericalCoor3D(qs.r,SphericalCoor3D(*vi).phi,SphericalCoor3D(*vi).theta) );
		}
	}	
	else if (m_vectors=="mcphiacos") {
		int num=0;
		while (num<m_resolution) {
			double p = 2*M_PI*(rand()*1.0/RAND_MAX);
			double t = acos(2*(rand()*1.0/RAND_MAX)-1);		
			
			// take length from scattering vector and angles from unisphere			
			m_qvectors.push_back( SphericalCoor3D(qs.r,p,t) );

			num++;
		}
	}	
	else if (m_vectors=="mcdoublesqrt") {
		double x1,x2;
		while (num<m_resolution) {
			x1 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x2 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			double xl = powf(x1,2) + powf(x2,2);
			if ( xl >= 1.0 ) continue;

			double x = 2*x1* sqrt(1-xl);
			double y = 2*x2* sqrt(1-xl) / xl;
			double z = 1-2*xl;

			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			// take length from scattering vector and angles from unisphere			
			m_qvectors.push_back( SphericalCoor3D(qs.r,s.phi,s.theta) );

			num++;
		}
	}
	else if (m_vectors=="mcquaternion") {
		double x0,x1,x2,x3;
		
		while (num<m_resolution) {
			x0 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x1 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x2 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			x3 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
			double xl = powf(x0,2) + powf(x1,2) + powf(x2,2) + powf(x3,2);
			if ( xl >= 1.0 ) continue;

			double x = 2* (x1*x3+x0*x2) / xl;
			double y = 2* (x2*x3-x0*x1) / xl;
			double z = 2* (powf(x0,2)+powf(x3,2)-powf(x1,2)+powf(x2,2)) / xl;

			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			// take length from scattering vector and angles from unisphere			
			m_qvectors.push_back( SphericalCoor3D(qs.r,s.phi,s.theta) );
			num++;
		}
	}
	else if (m_vectors=="mcboostunisphere") {
		
		boost::mt19937 rng; // that's my random number generator
		rng.seed(m_seed);		
		boost::uniform_on_sphere<double> s(3); // that's my distribution
		boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<double> > mysphere(rng,s);

		while (num<m_resolution) {

			vector<double> r = mysphere();

			double x = r[0];
			double y = r[1];
			double z = r[2];
			
			SphericalCoor3D s = CartesianCoor3D(x,y,z);
			// take length from scattering vector and angles from unisphere			
			m_qvectors.push_back( SphericalCoor3D(qs.r,s.phi,s.theta) ); 
			
			num++;
		}
	}
	else {
		Err::Inst()->write("didn't recognize vectors-keyword: ");
		throw;
	}	
}


CylinderVectorUnfold::CylinderVectorUnfold(CartesianCoor3D q) {
	m_q = q;
}

void CylinderVectorUnfold::set_resolution(size_t resolution) {
	m_resolution = resolution;
}

void CylinderVectorUnfold::set_axis(CartesianCoor3D axis) {
	m_axis = axis;
}

void CylinderVectorUnfold::set_seed(uint32_t seed) {
	m_seed = seed;
}

void CylinderVectorUnfold::set_vectors(string vectors) {
	m_vectors = vectors;
}

vector<CartesianCoor3D>& CylinderVectorUnfold::qvectors() {
	return m_qvectors;
}

void CylinderVectorUnfold::execute() {
	
	srand(m_seed);
	
	const double M_2PI = 2*M_PI;
	const double radincr = (M_2PI*m_resolution)/360;			
	
	CartesianCoor3D o = m_axis;
	// make sure o is normalized;
	o = o / o.length();
	
	// get the part of the scattering vector perpenticular to the o- orientation
	CartesianCoor3D qparallel = (o*m_q)*o; 
	CartesianCoor3D qperpenticular = m_q - qparallel; 

	// resolution is in degree, coordinates are in rad
	
	double qperpenticular_l = qperpenticular.length();
	double qparallel_l = qparallel.length();
	
	// we can't unfold a line:
	if (qperpenticular_l==0.0)  { 
		m_qvectors.push_back( qparallel );		
		return;
	}
	else {
		
		// rotation about an arbitrary axis u
		// R = P + (I-P)*cos(phi)+Q*sin(phi)
		// w/ P= [[ ux*ux, ux*uy , ux*uz ], [ ux*uy, uy*uy, uy*uz], [ux*uz, uy*uz, uz*uz] ]
		// w/ Q = [ [0, -uz, uy], [uz, 0, -ux], [-uy, ux, 0] ]
		// Rx = ux*ux + (1-ux*ux)*cos(phi) , ux*uy* (1- cos(phi)) -uz*sin(phi), ux*uz - ux*uz*cos(phi) +uy*sin(phi)

		// here: u = o
		// q' = R * q
		
		if (m_vectors=="rasterlinear") {
			for (double phi=0;phi<M_2PI;phi+=radincr) {
				CartesianCoor3D qv = m_q;
				CartesianCoor3D qv2;	
				qv2.x = (o.x*o.x+(1-o.x*o.x)*cos(phi))*qv.x      + (o.x*o.y*(1-cos(phi))-o.z*sin(phi))*qv.y + (o.x*o.z*(1-cos(phi))+o.y*sin(phi))*qv.z;
				qv2.y = (o.x*o.y*(1-cos(phi))+o.z*sin(phi))*qv.x + (o.y*o.y+(1-o.y*o.y)*cos(phi))*qv.y      + (o.y*o.z*(1-cos(phi))-o.x*sin(phi))*qv.z;
				qv2.z = (o.x*o.z*(1-cos(phi))-o.y*sin(phi))*qv.x + (o.y*o.z*(1-cos(phi))+o.x*sin(phi))*qv.y + (o.z*o.z+(1-o.z*o.z)*cos(phi))*qv.z;

				// take length from scattering vector and angles from unisphere										
				m_qvectors.push_back( qv2 );	
			}	
		}
		else {
			Err::Inst()->write("vectors generate method not understood");
			throw;
		}

	}
}
