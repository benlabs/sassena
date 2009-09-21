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

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>

#include "log.hpp"
#include "parameters.hpp"
#include "timer.hpp"

using namespace std;

PerformanceAnalyzer::PerformanceAnalyzer(boost::mpi::communicator thisworld, Timer timer) {
	// aggregate all timer on head node
	boost::mpi::gather(thisworld,timer,m_alltimer,0);
	
	analyze();
}


void PerformanceAnalyzer::analyze() {
	m_supertimer.clear();
	for(size_t i = 0; i < m_alltimer.size(); ++i)
	{
		m_supertimer += m_alltimer[i];
	}
}

void PerformanceAnalyzer::report() {

	Timer& timer = m_supertimer;

	using boost::io::group;
		vector<string> keys = timer.keys();
		string mthead = (boost::format("%|17| |%|9| |%|12| |%|9| |%|9|") % "measure" % "total" % "count" % "mean" % "stddev").str();
		string mmhead = (boost::format("%|17| |%|9| |%|9|") % "measure" % "min" % "max").str();

		Info::Inst()->write("                                                                 ");
		Info::Inst()->write("                    Performance Analysis                         ");
		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(" mean and total runtimes:                                        ");				
		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(mthead);
		Info::Inst()->write("-----------------------------------------------------------------");
		
		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			boost::format form(" %|17|| %|9|| %|12|| %|9|| %|9|");
			form % *ki;
			form % group(setprecision(3),scientific, timer.sum(*ki)            );
			form % timer.count(*ki) ;
			form % group(setprecision(3),scientific, timer.mean(*ki)           );
			form % group(setprecision(3),scientific, sqrt(timer.variance(*ki)) );
			Info::Inst()->write(form.str());
		}
		Info::Inst()->write("-----------------------------------------------------------------");

		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(" watermarks:                                                     ");				
		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(mmhead);
		Info::Inst()->write("-----------------------------------------------------------------");

		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
			boost::format form(" %|17|| %|9|| %|9|");
			form % *ki;
			form % group(setprecision(3) ,scientific, timer.min(*ki) );
			form % group(setprecision(3) ,scientific, timer.max(*ki) );
			Info::Inst()->write(form.str());
		}
		Info::Inst()->write("-----------------------------------------------------------------");

}