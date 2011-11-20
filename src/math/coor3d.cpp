/** \file
This file contains definitions for coordinate vector types for different coordinate systems.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "math/coor3d.hpp"

// standard header
#include <cmath>
#include <iostream>

// special library headers

// other headers

using namespace std;

/* CARTESIAN COOR3D CLASS */
bool CartesianCoor3D::operator<(const CartesianCoor3D& that) const{
	if (this->x < that.x)
		return true;

	if (this->x > that.x)
		return false;

	if (this->y < that.y)
		return true;

	if (this->y > that.y)
		return false;

	if (this->z < that.z)
		return true;

	if (this->z > that.z)
		return false;

	return false;
}
//conversion constructor cylinder -> cartesian
CartesianCoor3D::CartesianCoor3D(CylinderCoor3D cc) {
	x = cc.r*cos(cc.phi);
	y = cc.r*sin(cc.phi);
	z = cc.z;
}

//conversion constructor spherical -> cartesian
CartesianCoor3D::CartesianCoor3D(SphericalCoor3D cc) {
	x = cc.r*sin(cc.theta)*cos(cc.phi);
	y = cc.r*sin(cc.theta)*sin(cc.phi);
	z = cc.r*cos(cc.theta);
}

coor2_t CartesianCoor3D::length() { 
	return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}


ostream& operator<<(ostream& os, const CartesianCoor3D& cc) { 
	return os << "(x=" << cc.x << ",y=" << cc.y << ",z=" << cc.z << ")";
}


CartesianCoor3D& CartesianCoor3D::operator=(const CartesianCoor3D& that) {
	if (this != &that) { x=that.x; y=that.y; z=that.z; }
	return *this;
}

CartesianCoor3D CartesianCoor3D::operator-(const CartesianCoor3D& that) {
	return CartesianCoor3D(x-that.x,y-that.y,z-that.z);
}

CartesianCoor3D CartesianCoor3D::operator+(const CartesianCoor3D& that) {
	return CartesianCoor3D(x+that.x,y+that.y,z+that.z);
}

coor2_t CartesianCoor3D::operator*(const CartesianCoor3D& that) {
	return (x*that.x+y*that.y+z*that.z);
}

CartesianCoor3D CartesianCoor3D::cross_product(const CartesianCoor3D& that) {
	return CartesianCoor3D(y*that.z-z*that.y,z*that.x-x*that.z,x*that.y-y*that.x);
}

CartesianCoor3D operator*(const coor2_t lambda, const CartesianCoor3D& that) {
	return CartesianCoor3D(lambda*that.x,lambda*that.y,lambda*that.z);
}

CartesianCoor3D operator*(const CartesianCoor3D& that,const coor2_t lambda) {
	return CartesianCoor3D(lambda*that.x,lambda*that.y,lambda*that.z);
}

CartesianCoor3D operator/(const CartesianCoor3D& that, const coor2_t lambda) {
	return CartesianCoor3D(that.x/lambda,that.y/lambda,that.z/lambda);
}




/* CYLINDER COOR3D CLASS */
//conversion constructor cartesian -> cylinder
CylinderCoor3D::CylinderCoor3D(coor2_t v1,coor2_t v2,coor2_t v3) {

	if (v1<0) throw; // that's illegal

	r   = v1;
	phi = v2;
	z   = v3;
}

//conversion constructor cartesian -> cylinder
CylinderCoor3D::CylinderCoor3D(CartesianCoor3D cc) {

	r = sqrt(pow(cc.x,2)+pow(cc.y,2));
	
	if (cc.x!=0.0) {
		phi = atan(cc.y/cc.x);
		if (cc.x<0.0) { phi = sign(M_PI,cc.y)+phi; }
	}
	else if (cc.y!=0.0) {
		phi = sign(M_PI_2,cc.y);
	}
	else {
		phi = 0.0;
	}
	
	phi = phi<0 ? 2*M_PI + phi : phi;
	
	if ( (phi<0) || (phi>2*M_PI)) {
		cerr << "PHI OUT OF BOUND: " << r << ", " << phi << endl;
		cerr << cc << endl;
		
		throw;
	}
	z = cc.z;
}


//conversion constructor spherical -> cylinder
CylinderCoor3D::CylinderCoor3D(SphericalCoor3D cc) {
	//fixme: costly
	CartesianCoor3D c1(cc);
	CylinderCoor3D c2(c1);
	r = c2.r; phi = c2.phi; z = c2.z;
}


ostream& operator<<(ostream& os, const CylinderCoor3D& cc) { 
	return os << "(r=" << cc.r << ",phi=" << cc.phi << ",z=" << cc.z << ")";
}

CylinderCoor3D& CylinderCoor3D::operator=(const CylinderCoor3D& that) {
	if (this != &that) { r=that.r; phi=that.phi; z=that.z; }
	return *this;
}

CylinderCoor3D CylinderCoor3D::operator-(const CylinderCoor3D& that) {
	return CylinderCoor3D(r-that.r,phi-that.phi,z-that.z);
}

/* SPHERICAL COOR3D CLASS */

//conversion constructor cartesian -> spherical
SphericalCoor3D::SphericalCoor3D(CartesianCoor3D cc) {	
	r = sqrt(pow(cc.x,2)+pow(cc.y,2)+pow(cc.z,2));
    
    if (r==0) {
        theta = 0;
        phi=0;
        return;
    }

	theta = acos(cc.z/r);
	if (theta<0) throw;
	if (theta>M_PI) throw;	
//	theta = theta<0 ? 2*M_PI + theta : theta;

//	if (cc.z!=0.0) {
//		if (cc.z<0.0) { theta += M_PI; }
//	}
//	else {
//		theta = M_PI_2;
//	}



	if (cc.x!=0.0) {
		phi = atan(cc.y/cc.x);
		if (cc.x<0.0) phi += M_PI;
		else if (cc.y<0.0) phi += 2*M_PI;
	}
	else if (cc.y!=0.0) {
		if (cc.y>0) phi=M_PI_2;
		if (cc.y<0) phi=3*M_PI_2;
	}
	else {
		phi = 0.0;
	}
	if (phi<0) throw;
	if (phi>2*M_PI) throw;	
//	phi = phi<0 ? 2*M_PI + phi : phi;

//  if (cc.x!=0.0) {
//  	phi = atan(cc.y/cc.x);
//		if (cc.x<0.0) { phi = sign(M_PI,cc.y)+phi; }
//  	if (cc.y<0.0) { phi += M_PI; }
//  }
//  else {
//  	phi = sign(M_PI_2,cc.y);
//	}
}

//conversion constructor cylinder -> spherical
SphericalCoor3D::SphericalCoor3D(CylinderCoor3D cc) {
	//fixme: costly
	CartesianCoor3D c1(cc);
	SphericalCoor3D c2(c1);
	r = c2.r; phi = c2.phi; theta = c2.theta;
}

ostream& operator<<(ostream& os, const SphericalCoor3D& cc) { 
	return os << "(r=" << cc.r << ",phi=" << cc.phi << ",theta=" << cc.theta << ")";
}

SphericalCoor3D& SphericalCoor3D::operator=(const SphericalCoor3D& that) {
	if (this != &that) { r=that.r; phi=that.phi; theta=that.theta;}
	return *this;
}

SphericalCoor3D SphericalCoor3D::operator-(const SphericalCoor3D& that) {
	return SphericalCoor3D(r-that.r,phi-that.phi,theta-that.theta);
}


// transformations

CartesianCoor3D rotate(CartesianCoor3D c,string axis,coor2_t rad) {
	CartesianCoor3D r;
	if (axis=="x") {
		double m[3][3];
		m[0][0] = 1; m[0][1] = 0; m[0][2] = 0;
		m[1][0] = 0; m[1][1] = cos(rad); m[1][2] = -sin(rad);
		m[2][0] = 0; m[2][1] = sin(rad); m[2][2] = cos(rad);
		
		r.x = m[0][0]*c.x + m[0][1]*c.y + m[0][2]*c.z;
		r.y = m[1][0]*c.x + m[1][1]*c.y + m[1][2]*c.z;
		r.z = m[2][0]*c.x + m[2][1]*c.y + m[2][2]*c.z;
	}
	if (axis=="y") {
		double m[3][3];
		m[0][0] = cos(rad); m[0][1] = 0; m[0][2] = sin(rad);
		m[1][0] = 0; m[1][1] = 1; m[1][2] = 0;
		m[2][0] = -sin(rad); m[2][1] = 0; m[2][2] = cos(rad);
		
		r.x = m[0][0]*c.x + m[0][1]*c.y + m[0][2]*c.z;
		r.y = m[1][0]*c.x + m[1][1]*c.y + m[1][2]*c.z;
		r.z = m[2][0]*c.x + m[2][1]*c.y + m[2][2]*c.z;		
		
	}
	if (axis=="z") {
		double m[3][3];
		m[0][0] = cos(rad); m[0][1] = -sin(rad); m[0][2] = 0;
		m[1][0] = sin(rad); m[1][1] = cos(rad); m[1][2] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 1;
		
		r.x = m[0][0]*c.x + m[0][1]*c.y + m[0][2]*c.z;
		r.y = m[1][0]*c.x + m[1][1]*c.y + m[1][2]*c.z;
		r.z = m[2][0]*c.x + m[2][1]*c.y + m[2][2]*c.z;		
		
	}
	return r;
}

CartesianVectorBase::CartesianVectorBase(CartesianCoor3D axis) {
    CartesianCoor3D ek(0,0,1);
    CartesianCoor3D ej(0,1,0);
    
    
    CartesianCoor3D ez = axis/axis.length();
    CartesianCoor3D ekez = ek.cross_product(ez);
    if (ekez.length()==0) {
        ekez = ej.cross_product(ez);
    }
    CartesianCoor3D er = ekez/ekez.length();
    CartesianCoor3D ezer = ez.cross_product(er);
    CartesianCoor3D ephi = ezer/ezer.length();
    
    base_.push_back(er);
    base_.push_back(ephi);
    base_.push_back(ez);
}

CartesianCoor3D& CartesianVectorBase::operator[](size_t index) {
    return base_.at(index);
}

// project the vector vec onto the base and return the 3 components along the three basis vectors
CartesianCoor3D CartesianVectorBase::project(CartesianCoor3D vec) {
    return CartesianCoor3D(vec*base_[0],vec*base_[1],vec*base_[2]);
}

// end of file
