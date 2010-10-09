/*
 *  coordinateset.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
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

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

enum CoordinateRepresentation { 
    CARTESIAN=10, SPHERICAL=20, CYLINDRICAL=30 
};

// This class is used by Frame to store selection specific coordinate arrays
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
	CoordinateSet(CoordinateSet& cs,Atomselection& cs_selection, Atomselection& sub_selection); 
	
	std::vector<coor2_t> c1; // x-coordinates
	std::vector<coor2_t> c2; // y-coordinates
	std::vector<coor2_t> c3; // z-coordinates
		
	size_t size() { return m_size; }
	
    CoordinateRepresentation get_representation() {return m_representation;}
};

// This class is used by Frame to store selection specific coordinate arrays
class CartesianCoordinateSet : public CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<CoordinateSet,CartesianCoordinateSet>(*this);
    }

public:
	CartesianCoordinateSet();
	CartesianCoordinateSet(CartesianCoordinateSet& cs,Atomselection& cs_selection, Atomselection& sub_selection); 
	CartesianCoordinateSet(Frame& frame,Atomselection& selection); 

    // implement transformations in Cartesian Space
	void translate(CartesianCoor3D trans);
	void translate(CartesianCoor3D trans, Atomselection& cs_selection, Atomselection& sub_selection);
	
    void rotate(CartesianCoor3D axis1,CartesianCoor3D axis2);
};

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

class CylindricalCoordinateSet : public CoordinateSet {
	friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<CoordinateSet,CylindricalCoordinateSet>(*this);
    }

public:
	CylindricalCoordinateSet();
	CylindricalCoordinateSet(CartesianCoordinateSet& cs); 
};

#endif

// end of file
