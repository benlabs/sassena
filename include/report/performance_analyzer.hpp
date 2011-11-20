/** \file
This file contains an performance analysis class, which takes several timers and prepares a formatted output.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include <boost/thread.hpp>

// other headers
#include "report/timer.hpp"

using namespace boost::accumulators;

/** 
Data entry used for performance analysis and contains timer information.
*/
struct PerformanceMeasure {
	double sum;
	double count;
	double mean;
	double variance;
	double min;
	double max;
};

/** 
Analyzes the information from a set of timers provided at construction and prepares a formatted output of the corresponding statistical information.
*/
class PerformanceAnalyzer {

//	Timer m_supertimer; // contains all timing information from all nodes
//	std::vector<Timer> m_alltimer;
	std::map<std::string, PerformanceMeasure > m_measures;
	
	void analyze();
	
public:
	PerformanceAnalyzer(boost::mpi::communicator thisworld,std::map<boost::thread::id,Timer>& timer);
	
	void report();	
    void report_relative(double totaltime);
};


#endif

// end of file 
