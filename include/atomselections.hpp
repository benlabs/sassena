/*
 *  atomselections.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef ATOMSELECTIONS_HPP_
#define ATOMSELECTIONS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <map>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>

// other headers
#include "atomselection.hpp"

class Atomselections : public std::map<std::string, Atomselection> {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::map<std::string,Atomselection> >(*this);
    }
	///////////////////
};

#endif
