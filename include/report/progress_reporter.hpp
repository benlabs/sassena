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
	
	double m_maximum;
	double m_minimum;
	double m_current;
	
	double m_nextreport;
	
	std::string m_timerlabel;

	double compute_nextreport();

public:
	ProgressReporter(double minimum,double maximum);
	
	void set(double current);
	void report();

};

#endif

// end of file
