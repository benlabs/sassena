/*
 *  scatterdata.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTERDATA_HPP_
#define SCATTERDATA_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

// other headers
#include "coor3d.hpp"

class ScatterdataKey : public std::pair<CartesianCoor3D,int> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::pair<CartesianCoor3D,int> >(*this);
    }
	///////////////////
	
public:
	ScatterdataKey() {}
	ScatterdataKey(CartesianCoor3D q, int frame) { this->first = q; this->second = frame; }	
	
	bool operator <(const ScatterdataKey that) const;
};

template <class V> class Scatterdata : public std::map<ScatterdataKey,V> {
	
	std::map<ScatterdataKey,std::pair<V,V> > average();

public:
	void write_plain(std::string fname,std::string format);
	void write_average(std::string fname,std::string format);
	void dump(std::streambuf* pof);
	
};

#endif
