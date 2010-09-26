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
#include "log.hpp"

using namespace std;

ProgressReporter::ProgressReporter(double minimum, double maximum) {
	// use internal timer to track elapsed time
	m_maximum = maximum;
	m_minimum = minimum;
	m_current = minimum;
	
	m_nextreport = minimum+0.001;
	
	m_timerlabel = "progress";
	m_timer.start(m_timerlabel);
}

void ProgressReporter::set(double c) {
	m_current = c;
	m_timer.stop(m_timerlabel);
	m_timer.start(m_timerlabel);
}

void ProgressReporter::report() {
	if (m_current >= m_nextreport) {
		m_nextreport = compute_nextreport();
		
		double fraction_done = (m_current-m_minimum) / (m_maximum-m_minimum);
		//average time since start: 
		double totaltime = m_timer.sum(m_timerlabel);
		
		double etotal = ( 1.0/fraction_done) * totaltime;
		double eta = etotal - totaltime;
	
        stringstream ifraction_done_strstr;
        ifraction_done_strstr.precision(3);
        ifraction_done_strstr.setf(std::ios::fixed);
        std::cout << "m_minimum" << m_minimum << " m_maximum" << m_maximum << " m_current" << m_current << std::endl;
        ifraction_done_strstr << 100*fraction_done;
        
		Info::Inst()->write(string("Progress: ") + ifraction_done_strstr.str() +string("%, ETOTAL: ")+to_s(etotal) +string(" s, ETA: ")+to_s(eta)+ string(" s"));
	}
}

double ProgressReporter::compute_nextreport() {
	
	// setup a delay pattern for the report: 1,2,3,5,10,20,30,40,50,60,70,80,90,100
	// double fraction = (m_current-m_minimum)*100.0 / (m_maximum-m_minimum);

    vector<double> delays;
    delays.push_back(0.001);
    delays.push_back(0.01);
    delays.push_back(0.1);
    delays.push_back(1);
    delays.push_back(2);
    delays.push_back(3);
    delays.push_back(5);
    delays.push_back(10);
    delays.push_back(20);
    delays.push_back(30);
    delays.push_back(40);
    delays.push_back(50);
    delays.push_back(60);
    delays.push_back(70);
    delays.push_back(80);
    delays.push_back(90);
    delays.push_back(100);

	vector<double> next_reports;

	for(size_t i = 0; i < delays.size(); ++i)
	{
		double threshold = (m_maximum - m_minimum)*delays[i]/100.0;
		double nreport = threshold + m_minimum;
		next_reports.push_back(nreport);
	}

    double nextreport = m_maximum;
    
	for(size_t i = 0; i < next_reports.size(); ++i)
	{
		if (m_current < next_reports[i]) {
		    nextreport = next_reports[i];
            break;
	    }
	}
    cout << "m_current=" << m_current << ",nextrep=" << nextreport << endl;
	
	return nextreport; // default, if m_current is not within next_reports
}

// end of file
