/*
 *  coor3d.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef COOR3D_HPP_
#define COOR3D_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>

// special library headers
#include <boost/serialization/access.hpp>

// other headers

class CartesianCoor3D;
class CylinderCoor3D ;
class SphericalCoor3D;

typedef std::pair< CartesianCoor3D,CartesianCoor3D> cartrect;

//helper funcions:
inline double sign(double a,double b) { return (b<0.0) ? -a : a; }
inline float sign(float a,float b) { return (b<0.0) ? -a : a; }

class CartesianCoor3D {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & x;
		ar & y;
		ar & z;		
    }
	///////////////////
public:
	double x,y,z;
	
	CartesianCoor3D() : x(0), y(0), z(0) {}	
	CartesianCoor3D(double v1,double v2,double v3) : x(v1), y(v2), z(v3) {}
	// conversion constructors:
	CartesianCoor3D(const CylinderCoor3D cc);
	CartesianCoor3D(const SphericalCoor3D cc);
	
	double length();
	
	friend std::ostream& operator<<(std::ostream& os, const CartesianCoor3D& cc);
	
	CartesianCoor3D& operator=(const CartesianCoor3D& that);
	CartesianCoor3D operator-(const CartesianCoor3D& that);
	CartesianCoor3D operator+(const CartesianCoor3D& that);	
	CartesianCoor3D cross_product(const CartesianCoor3D& that);
	double operator*(const CartesianCoor3D& that);
	
	// for use in maps only!
	bool operator<(const CartesianCoor3D& that) const;
		
	~CartesianCoor3D() {}
};

CartesianCoor3D operator*(const double lambda, const CartesianCoor3D& that);
CartesianCoor3D operator*(const CartesianCoor3D& that,const double lambda);
CartesianCoor3D operator/(const CartesianCoor3D& that,const double lambda);

// cylinder coords have a range:
// 0 <= r
// 0 <= phi < 2 M_PI
class CylinderCoor3D {
public:
	double r,phi,z;
	
	CylinderCoor3D() : r(0), phi(0), z(0) {}	
	CylinderCoor3D(double v1,double v2,double v3);
	// conversion constructors:
	CylinderCoor3D(const CartesianCoor3D cc);
	CylinderCoor3D(const SphericalCoor3D cc);	
	
	friend std::ostream& operator<<(std::ostream& os, const CylinderCoor3D& cc);
	
	CylinderCoor3D& operator=(const CylinderCoor3D& that);
	CylinderCoor3D operator-(const CylinderCoor3D& that);
	
	~CylinderCoor3D() {}	
};


class SphericalCoor3D {
public:	
	double r,phi,theta;
	
	SphericalCoor3D() : r(0), phi(0), theta(0) {}	
	SphericalCoor3D(double v1,double v2,double v3) : r(v1), phi(v2), theta(v3) {}
	// conversion constructors:
	SphericalCoor3D(const CartesianCoor3D cc);
	SphericalCoor3D(const CylinderCoor3D cc);
	
	friend std::ostream& operator<<(std::ostream& os, const SphericalCoor3D& cc);
		
	SphericalCoor3D& operator=(const SphericalCoor3D& that);
	SphericalCoor3D operator-(const SphericalCoor3D& that);

	~SphericalCoor3D() {}
};

CartesianCoor3D rotate(CartesianCoor3D,std::string axis,double rad);

#endif
