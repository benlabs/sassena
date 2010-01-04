/*
 *  atomselection.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SAMPLE__ATOMSELECTION_HPP_
#define SAMPLE__ATOMSELECTION_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>

class Atomselection {
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & indexes;
    }
	///////////////////
public:	
    std::vector<size_t> indexes;    
    
	Atomselection& operator+=(const Atomselection& );
	Atomselection& operator-=(const Atomselection& );
	
	void add(size_t index);
	void remove(size_t index);
	
	std::vector<Atomselection> slice(size_t parts);
};

#endif

// end of file
