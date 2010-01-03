/*
 *  progress_reporter.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "report/progress_reporter.hpp"
#include "control.hpp"

using namespace std;

ProgressReporter::ProgressReporter(size_t maximum, size_t minimum) {
	// use internal timer to track elapsed time
	m_maximum = maximum;
	m_minimum = minimum;
	m_current = minimum;
	
	m_nextreport = minimum+1;
	
	m_timerlabel = "progress";
	m_timer.start(m_timerlabel);
}

void ProgressReporter::set(size_t c) {
	m_current = c;
	m_timer.stop(m_timerlabel);
	m_timer.start(m_timerlabel);
}

void ProgressReporter::report() {
	if (m_current >= m_nextreport) {
		m_nextreport = compute_nextreport();
		
		double fraction_done = (m_current-m_minimum)*1.0 / (m_maximum-m_minimum);
		//average time since start: 
		double totaltime = m_timer.sum(m_timerlabel);
		
		double etotal = ( 1.0/fraction_done) * totaltime;
		double eta = etotal - totaltime;
	
		size_t ifraction_done = long(10000*fraction_done)/100;

		Info::Inst()->write(string("Progress: ") + to_s(ifraction_done) +string("%, ETOTAL: ")+to_s(etotal) +string(" s, ETA: ")+to_s(eta)+ string(" s"));
	}
}

size_t ProgressReporter::compute_nextreport() {
	
	// setup a delay pattern for the report: 1,2,3,5,10,20,30,40,50,60,70,80,90,100
	double fraction = (m_current-m_minimum)*100.0 / (m_maximum-m_minimum);

	double delaysa[] = {1,2,3,5,10,20,30,40,50,60,70,80,90,100};
	vector<double> delays(delaysa,delaysa+13);

	vector<size_t> next_reports;

	for(size_t i = 0; i < delays.size(); ++i)
	{
		double threshold = (m_maximum - m_minimum)*delays[i]/100.0;
		size_t nreport = threshold + m_minimum;
		next_reports.push_back(nreport);
	}

	for(size_t i = 0; i < next_reports.size(); ++i)
	{
		if (m_current < next_reports[i]) return next_reports[i];
	}
	
	return m_maximum; // default
}
