/*
 *  timer.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "sample.hpp"

using namespace boost::accumulators;

class Timer {
	typedef boost::accumulators::accumulator_set<double,features<tag::min,tag::max,tag::mean,tag::variance,tag::sum,tag::count> > times_type;
	std::map<std::string,times_type > times;
	
	std::map<std::string,timeval> starttimes;
	std::map<std::string,bool> states;
	
	double t_diff(timeval start, timeval end);
	

public:			
	std::vector<std::string> keys();
	
	void start(std::string tk);
	void stop(std::string tk);

	double sum(std::string tk);
	double min(std::string tk);
	double max(std::string tk);	
	double mean(std::string tk);
	double variance(std::string tk);	
	double count(std::string tk);		
};

#endif TIMER_HPP_
