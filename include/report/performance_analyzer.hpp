/*
 *  performance_analyzer.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef REPORT__PERFORMANCE_ANALYZER_HPP_
#define REPORT__PERFORMANCE_ANALYZER_HPP_

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
#include <boost/mpi.hpp>

// other headers
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "sample/sample.hpp"
#include "report/timer.hpp"

using namespace boost::accumulators;

struct PerformanceMeasure {
	double sum;
	double count;
	double mean;
	double variance;
	double min;
	double max;
};

class PerformanceAnalyzer {

//	Timer m_supertimer; // contains all timing information from all nodes
//	std::vector<Timer> m_alltimer;
	std::map<std::string, PerformanceMeasure > m_measures;
	
	void analyze();
	
public:
	PerformanceAnalyzer(boost::mpi::communicator thisworld, Timer timer);
	
	void report();	
    void report_relative(double totaltime);
};


#endif

// end of file 
