/** \file
This file contains a class which defines coordinates based on a coordinate system.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef SAMPLE__COORDINATESET_HPP_
#define SAMPLE__COORDINATESET_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "sample/frame.hpp"

//forward declaration...

/** 
Defines the three possible coordinate representations
*/
enum CoordinateRepresentation { 
    CARTESIAN=10, SPHERICAL=20, CYLINDRICAL=30 
};

/** 
A set of coordinates with an associated representation
*/
class CoordinateSet {
protected:
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {        
		ar & c1;
		ar & c2;
		ar & c3;
		ar & m_size;
        ar & m_representation;
    }
	size_t m_size;
    CoordinateRepresentation m_representation;
public:
    CoordinateSet();
	CoordinateSet(CoordinateSet& cs,IAtomselection* pcs_selection, IAtomselection* psub_selection); 
	
	std::vector<coor2_t> c1; // x-coordinates
	std::vector<coor2_t> c2; // y-coordinates
	std::vector<coor2_t> c3; // z-coordinates
		
	size_t size() { return m_size; }
	
    CoordinateRepresentation get_representation() {return m_representation;}
};

/** 
Specialized coordinate set which provides addtional functions for translation and rotation.
*/
class CartesianCoordinateSet : public CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<CoordinateSet,CartesianCoordinateSet>(*this);
    }

public:
	CartesianCoordinateSet();
	CartesianCoordinateSet(CartesianCoordinateSet& cs,IAtomselection* pcs_selection, IAtomselection* psub_selection); 
	CartesianCoordinateSet(Frame& frame,IAtomselection* selection); 

    // implement transformations in Cartesian Space
	void translate(CartesianCoor3D trans);
	void translate(CartesianCoor3D trans, IAtomselection* pcs_selection, IAtomselection* psub_selection);
	void transform(boost::numeric::ublas::matrix<double> T, IAtomselection* pcs_selection, IAtomselection* psub_selection);
	
    void rotate(CartesianCoor3D axis1,CartesianCoor3D axis2);
};

/** 
Specialized coordinate set for spherical representation
*/
class SphericalCoordinateSet : public CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<CoordinateSet,SphericalCoordinateSet>(*this);
    }

public:
	SphericalCoordinateSet();
	SphericalCoordinateSet(CartesianCoordinateSet& cs); 
};

/** 
Specialized coordinate set for cylindrical representation
*/
class CylindricalCoordinateSet : public CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<CoordinateSet,CylindricalCoordinateSet>(*this);
        ar & axis_;
    }
    CartesianCoor3D axis_;
public:
	CylindricalCoordinateSet();
	CylindricalCoordinateSet(CartesianCoordinateSet& cs, CartesianCoor3D axis); 
	
    CartesianCoor3D get_axis() { return axis_;}
    void set_axis(CartesianCoor3D axis) {axis_ = axis;}
    
};

#endif

// end of file
