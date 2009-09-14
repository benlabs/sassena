/*
 *  coordinate_sets.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef COORDINATESETS_HPP_
#define COORDINATESETS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "atomselection.hpp"
#include "sample.hpp"
#include "frames.hpp"
#include "coordinate_set.hpp"

//forward declaration...

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

// This class is used by Scatterdevices to manage coordinate sets
class CoordinateSets {

	std::map<size_t,CoordinateSet*> setcache;

	Sample* p_sample;
	Atomselection* p_selection;
	
	size_t currentframe_i;

public:
	CoordinateSets();
	~CoordinateSets() ;

	CoordinateSet& current();
	CoordinateSet& load(size_t frame);	
	void clear_cache();
	
	// use these to initialize the coordinate set:
	void set_selection(Atomselection& selection);
	void set_sample(Sample& sample);
	
	Atomselection& get_selection();
};

#endif
// end of file