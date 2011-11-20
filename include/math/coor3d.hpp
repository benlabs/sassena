/** \file
This file contains definitions for coordinate vector types for different coordinate systems.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef MATH__COOR3D_HPP_
#define MATH__COOR3D_HPP_

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
inline float sign(float a,float b) { return (b<0.0) ? -a : a; }

/** 
Type class which represents coordinates in cartesian space. Allows transformation into other coordinate representations and implements some basic linear algebra.
*/
class CartesianCoor3D {
protected:
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
	coor2_t x,y,z;
	
	CartesianCoor3D() : x(0), y(0), z(0) {}	
	CartesianCoor3D(coor2_t v1,coor2_t v2,coor2_t v3) : x(v1), y(v2), z(v3) {}
	// conversion constructors:
	CartesianCoor3D(const CylinderCoor3D cc);
	CartesianCoor3D(const SphericalCoor3D cc);
	
	coor2_t length();
	
	friend std::ostream& operator<<(std::ostream& os, const CartesianCoor3D& cc);
	
	CartesianCoor3D& operator=(const CartesianCoor3D& that);
	CartesianCoor3D operator-(const CartesianCoor3D& that);
	CartesianCoor3D operator+(const CartesianCoor3D& that);	
	CartesianCoor3D cross_product(const CartesianCoor3D& that);
	coor2_t operator*(const CartesianCoor3D& that);
	
	// for use in maps only!
	bool operator<(const CartesianCoor3D& that) const;
		
	~CartesianCoor3D() {}
};

CartesianCoor3D operator*(const coor2_t lambda, const CartesianCoor3D& that);
CartesianCoor3D operator*(const CartesianCoor3D& that,const coor2_t lambda);
CartesianCoor3D operator/(const CartesianCoor3D& that,const coor2_t lambda);

/** 
Type class which represents coordinates in cylinder space. Allows transformation into other coordinate representations and implements some basic linear algebra.

Cylinder coords have a range:
- r >= 0
- 0 <= phi < 2 M_PI
*/
class CylinderCoor3D {
public:
	coor2_t r,phi,z;
	
	CylinderCoor3D() : r(0), phi(0), z(0) {}	
	CylinderCoor3D(coor2_t v1,coor2_t v2,coor2_t v3);
	// conversion constructors:
	CylinderCoor3D(const CartesianCoor3D cc);
	CylinderCoor3D(const SphericalCoor3D cc);	
	
	friend std::ostream& operator<<(std::ostream& os, const CylinderCoor3D& cc);
	
	CylinderCoor3D& operator=(const CylinderCoor3D& that);
	CylinderCoor3D operator-(const CylinderCoor3D& that);
	
	~CylinderCoor3D() {}	
};

/** 
Type class which represents coordinates in spherical space. Allows transformation into other coordinate representations and implements some basic linear algebra.

Spherical coords have a range:
- r >= 0
- 0 <= phi < 2 M_PI
- 0 <= theta < M_PI
*/
class SphericalCoor3D {
public:	
	coor2_t r,phi,theta;
	
	SphericalCoor3D() : r(0), phi(0), theta(0) {}	
	SphericalCoor3D(coor2_t v1,coor2_t v2,coor2_t v3) : r(v1), phi(v2), theta(v3) {}
	// conversion constructors:
	SphericalCoor3D(const CartesianCoor3D cc);
	SphericalCoor3D(const CylinderCoor3D cc);
	
	friend std::ostream& operator<<(std::ostream& os, const SphericalCoor3D& cc);
		
	SphericalCoor3D& operator=(const SphericalCoor3D& that);
	SphericalCoor3D operator-(const SphericalCoor3D& that);

	~SphericalCoor3D() {}
};

CartesianCoor3D rotate(CartesianCoor3D,std::string axis,coor2_t rad);

/** 
Type class which represents a vector base (3 orthonormal vectors) for cartesian coordinates. Can be constructed out of thin air or from partial vectors.
*/
class CartesianVectorBase {
    std::vector<CartesianCoor3D> base_;
public:
      CartesianVectorBase() {}
      CartesianVectorBase(CartesianCoor3D axis);
      std::vector<CartesianCoor3D> get_base() {return base_;}
      CartesianCoor3D& operator[](size_t index);
      CartesianCoor3D project(CartesianCoor3D vec);
};
#endif

// end of file
