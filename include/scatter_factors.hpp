/*
 *  scatter_factors.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_FACTORS_HPP_
#define SCATTER_FACTORS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "sample.hpp"

//forward declaration...

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

// This class is used by Scatterdevices to recalculate the q dependent scattering factors
class ScatterFactors {
	Sample* p_sample;
	Atomselection* p_selection;
	
	bool m_background;
	
	std::vector<double> factors;
public:
	ScatterFactors();
	
	// use these to initialize the coordinate set:
	void set_selection(Atomselection& selection);
	void set_sample(Sample& sample);
	
	Atomselection& get_selection();
	
	double get(size_t atomselectionindex);	
	std::vector<double>& get_all();
	
	void update(CartesianCoor3D q);
	void set_background(bool status);
};

#endif

// end of file
