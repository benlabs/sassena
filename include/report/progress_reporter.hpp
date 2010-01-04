/*
 *  progress_reporter.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef REPORT__PROGRESS_REPORTER_HPP_
#define REPORT__PROGRESS_REPORTER_HPP_

// common header
#include "common.hpp"

// standard header
#include <map>
#include <string>
#include <vector>

// special library headers

// other headers
#include "report/timer.hpp"

class ProgressReporter {

	Timer m_timer;
	
	size_t m_maximum;
	size_t m_minimum;
	size_t m_current;
	
	size_t m_nextreport;
	
	std::string m_timerlabel;

	size_t compute_nextreport();

public:
	ProgressReporter(size_t maximum,size_t minimum);
	
	void set(size_t current);
	void report();

};

#endif

// end of file
