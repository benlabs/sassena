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

#include "control.hpp"
#include "timer.hpp"

using namespace std;

PerformanceAnalyzer::PerformanceAnalyzer(boost::mpi::communicator thisworld, Timer timer) {
	// aggregate all timer on head node
//	boost::mpi::gather(thisworld,timer,m_alltimer,0);
	

	// we need: sum , count, mean , variance, min, max
	double my_sum, my_count, my_mean, my_variance, my_min, my_max;
	vector<double> all_sum, all_count, all_mean, all_variance, all_min, all_max;

	int my_keyflag;
	vector<int> all_keyflags;
	
	//	vector<double> all_sum, all_count, all_mean, all_variance, all_min, all_max;
	all_sum.resize(thisworld.size());
	all_count.resize(thisworld.size());
	all_mean.resize(thisworld.size());
	all_variance.resize(thisworld.size());
	all_min.resize(thisworld.size());
	all_max.resize(thisworld.size());

	all_keyflags.resize(thisworld.size());

	// first negotiate keys 
	vector<string> my_keys = timer.keys();
	vector< vector<string> > all_keys(thisworld.size());

	boost::mpi::gather(thisworld ,my_keys,all_keys,0);
	
	vector<string> common_keys;
	if (thisworld.rank()==0) {
		for(size_t i = 0; i < all_keys.size(); ++i)
		{	
			for(size_t j = 0; j < all_keys[i].size(); ++j)
			{
				if (find(common_keys.begin(),common_keys.end(),all_keys[i][j])==common_keys.end()) common_keys.push_back(all_keys[i][j]);
			}
		}		
	}
	
	
	boost::mpi::broadcast(thisworld,common_keys,0);
		
	for(vector<string>::iterator ti = common_keys.begin(); ti != common_keys.end(); ++ti)
	{
		PerformanceMeasure& pm = m_measures[*ti];
		pm.sum = 0;
		pm.count = 0;
		pm.mean = 0;
		pm.variance = 0;
		pm.min = 0;
		pm.max = 0;

		if ( find(my_keys.begin(),my_keys.end(),*ti)==my_keys.end() ) my_keyflag = 0; else my_keyflag = 1;

		if (my_keyflag == 1) {
			my_sum = timer.sum(*ti);
			my_count = timer.count(*ti);
			my_mean = timer.mean(*ti);
			my_variance = timer.variance(*ti);
			my_min = timer.min(*ti);
			my_max = timer.max(*ti);			
		}
		
		boost::mpi::gather(thisworld ,my_keyflag  , all_keyflags ,0);		
		boost::mpi::gather(thisworld ,my_sum      , all_sum      ,0);
		boost::mpi::gather(thisworld ,my_count    , all_count    ,0);
		boost::mpi::gather(thisworld ,my_mean     , all_mean     ,0);
		boost::mpi::gather(thisworld ,my_variance , all_variance ,0);
		boost::mpi::gather(thisworld ,my_min      , all_min      ,0);
		boost::mpi::gather(thisworld ,my_max      , all_max      ,0);

		if (thisworld.rank()==0) {
			double weightedmean=0;
			double weightedvariance=0;
			double totalmin = 1e80;
			double totalmax = -1;

			for(size_t i = 0; i < thisworld.size(); ++i)
			{
				if ( all_keyflags[i] == 1 ) {
					pm.sum     += all_sum[i];
					pm.count   += all_count[i];
					weightedmean    += all_mean[i]*all_count[i];
					weightedvariance += all_variance[i]*all_count[i];
					if (totalmin>all_min[i]) totalmin = all_min[i];
					if (totalmax<all_max[i]) totalmax = all_max[i];
				}
			}
			pm.mean = weightedmean / pm.count;
			pm.variance = weightedvariance / pm.count;
			pm.min = totalmin;
			pm.max = totalmax;
		}
	}
	// analyze();
}


void PerformanceAnalyzer::analyze() {
//	m_supertimer.clear();
//	for(size_t i = 0; i < m_alltimer.size(); ++i)
//	{
//		m_supertimer += m_alltimer[i];
//	}
}


void PerformanceAnalyzer::report() {
	using boost::io::group;
	string mthead = (boost::format("%|17| |%|9| |%|12| |%|9| |%|9|") % "measure" % "total" % "count" % "mean" % "stddev").str();
	string mmhead = (boost::format("%|17| |%|9| |%|9|") % "measure" % "min" % "max").str();

	Info::Inst()->write("                                                                 ");
	Info::Inst()->write("                    Performance Analysis                         ");
	Info::Inst()->write("-----------------------------------------------------------------");
	Info::Inst()->write(" mean and total runtimes:                                        ");				
	Info::Inst()->write("-----------------------------------------------------------------");
	Info::Inst()->write(mthead);
	Info::Inst()->write("-----------------------------------------------------------------");
		
	for (map<string,PerformanceMeasure>::iterator mi=m_measures.begin();mi!=m_measures.end();mi++) {
			boost::format form(" %|17|| %|9|| %|12|| %|9|| %|9|");
			form % mi->first;
			form % group(setprecision(3),scientific, mi->second.sum            );
			form % mi->second.count ;
			form % group(setprecision(3),scientific, mi->second.mean           );
			form % group(setprecision(3),scientific, sqrt(mi->second.variance) );
			Info::Inst()->write(form.str());
		}
		Info::Inst()->write("-----------------------------------------------------------------");

		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(" watermarks:                                                     ");				
		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(mmhead);
		Info::Inst()->write("-----------------------------------------------------------------");

	for (map<string,PerformanceMeasure>::iterator mi=m_measures.begin();mi!=m_measures.end();mi++) {
		boost::format form(" %|17|| %|9|| %|9|");
		form % mi->first;
		form % group(setprecision(3) ,scientific, mi->second.min );
		form % group(setprecision(3) ,scientific, mi->second.max );
		Info::Inst()->write(form.str());
	}
	Info::Inst()->write("-----------------------------------------------------------------");

}


void PerformanceAnalyzer::report_relative(double totaltime) {
	using boost::io::group;
	string mthead = (boost::format("%|17| |%|9| |%|12| |%|9| |%|9|") % "measure" % "total" % "count" % "mean" % "stddev").str();
	string mmhead = (boost::format("%|17| |%|9| |%|9|") % "measure" % "min" % "max").str();

	Info::Inst()->write("                                                                 ");
	Info::Inst()->write("         Performance Analysis (relative to total time * nodes)   ");
	Info::Inst()->write("-----------------------------------------------------------------");
	Info::Inst()->write(" mean and total runtimes:                                        ");				
	Info::Inst()->write("-----------------------------------------------------------------");
	Info::Inst()->write(mthead);
	Info::Inst()->write("-----------------------------------------------------------------");
		
	for (map<string,PerformanceMeasure>::iterator mi=m_measures.begin();mi!=m_measures.end();mi++) {
			boost::format form(" %|17|| %|9|| %|12|| %|9|| %|9|");
			form % mi->first;
			form % group(setprecision(3),scientific, mi->second.sum / totaltime           );
			form % mi->second.count ;
			form % group(setprecision(3),scientific, mi->second.mean / totaltime          );
			form % group(setprecision(3),scientific, sqrt(mi->second.variance / totaltime) );
			Info::Inst()->write(form.str());
		}
		Info::Inst()->write("-----------------------------------------------------------------");

		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(" watermarks:                                                     ");				
		Info::Inst()->write("-----------------------------------------------------------------");
		Info::Inst()->write(mmhead);
		Info::Inst()->write("-----------------------------------------------------------------");

	for (map<string,PerformanceMeasure>::iterator mi=m_measures.begin();mi!=m_measures.end();mi++) {
		boost::format form(" %|17|| %|9|| %|9|");
		form % mi->first;
		form % group(setprecision(3) ,scientific, mi->second.min / totaltime);
		form % group(setprecision(3) ,scientific, mi->second.max / totaltime);
		Info::Inst()->write(form.str());
	}
	Info::Inst()->write("-----------------------------------------------------------------");

}

//void PerformanceAnalyzer::report2() {
//
////	Timer& timer = m_supertimer;
//
//	using boost::io::group;
////		vector<string> keys = timer.keys();
//		vector<string> keys = measures.keys();
//		string mthead = (boost::format("%|17| |%|9| |%|12| |%|9| |%|9|") % "measure" % "total" % "count" % "mean" % "stddev").str();
//		string mmhead = (boost::format("%|17| |%|9| |%|9|") % "measure" % "min" % "max").str();
//
//		Info::Inst()->write("                                                                 ");
//		Info::Inst()->write("                    Performance Analysis                         ");
//		Info::Inst()->write("-----------------------------------------------------------------");
//		Info::Inst()->write(" mean and total runtimes:                                        ");				
//		Info::Inst()->write("-----------------------------------------------------------------");
//		Info::Inst()->write(mthead);
//		Info::Inst()->write("-----------------------------------------------------------------");
//		
//		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
//			boost::format form(" %|17|| %|9|| %|12|| %|9|| %|9|");
//			form % *ki;
//			form % group(setprecision(3),scientific, timer.sum(*ki)            );
//			form % timer.count(*ki) ;
//			form % group(setprecision(3),scientific, timer.mean(*ki)           );
//			form % group(setprecision(3),scientific, sqrt(timer.variance(*ki)) );
//			Info::Inst()->write(form.str());
//		}
//		Info::Inst()->write("-----------------------------------------------------------------");
//
//		Info::Inst()->write("-----------------------------------------------------------------");
//		Info::Inst()->write(" watermarks:                                                     ");				
//		Info::Inst()->write("-----------------------------------------------------------------");
//		Info::Inst()->write(mmhead);
//		Info::Inst()->write("-----------------------------------------------------------------");
//
//		for (vector<string>::iterator ki=keys.begin();ki!=keys.end();ki++) {
//			boost::format form(" %|17|| %|9|| %|9|");
//			form % *ki;
//			form % group(setprecision(3) ,scientific, timer.min(*ki) );
//			form % group(setprecision(3) ,scientific, timer.max(*ki) );
//			Info::Inst()->write(form.str());
//		}
//		Info::Inst()->write("-----------------------------------------------------------------");
//
//}