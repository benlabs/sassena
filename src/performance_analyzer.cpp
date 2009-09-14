/*
 *  performance_analyzer.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "performance_analyzer.hpp"

#include <deque>

#include <iostream>
#include <iomanip>

#include "log.hpp"
#include "parameters.hpp"
#include "timer.hpp"

using namespace std;

PerformanceAnalyzer::PerformanceAnalyzer(boost::mpi::communicator thisworld, Timer timer) {
	// aggregate all timer on head node
	vector<Timer> alltimer;
	boost::mpi::gather(thisworld,timer,alltimer,0);
}

void PerformanceAnalyzer::report() {

	Timer& timer = m_supertimer;
	stringstream ss;
		vector<string> keys = timer.keys();
				
		ss << "INFO>> " << "                                                                         " << endl;
		ss << "INFO>> " << "                    Performance Analysis                                 " << endl;
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		ss << "INFO>> " << " mean and total runtimes:                                                " << endl;				
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		ss << "INFO>> ";
		ss << setw(31) << " measure |";
		ss << setw(12) << " total |";
		ss << setw(10) << " count |";
		ss << setw(10) << " mean |";
		ss << setw(10) << " stddev ";
		ss << endl;
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			ss << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			ss << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.sum(*ki)            << " |";	
			ss << setiosflags(ios::fixed) << setprecision(0) << setw(8) << timer.count(*ki)          << " |";		
			ss << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.mean(*ki)           << " |";
			ss << setiosflags(ios::fixed) << setprecision(3) << setw(8) << sqrt(timer.variance(*ki)) << endl;
		}
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;		

		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		ss << "INFO>> " << " watermarks:                                                " << endl;				
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;
		ss << "INFO>> ";
		ss << setw(31) << " measure |";
		ss << setw(12) << " min |";
		ss << setw(10) << " max ";
		ss << endl;
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			ss << "INFO>> " << setw(29) << *ki << " |" << "\t" ;
			ss << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.min(*ki)            << " |";	
			ss << setiosflags(ios::fixed) << setprecision(3) << setw(8) << timer.max(*ki)            << endl;
		}
		ss << "INFO>> " << "-------------------------------------------------------------------------" << endl;		
	
}