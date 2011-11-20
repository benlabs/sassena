/** \file
This file contains an performance analysis class, which takes several timers and prepares a formatted output.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "report/performance_analyzer.hpp"

#include <deque>

#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>

#include "control.hpp"
#include "log.hpp"
#include "mpi/wrapper.hpp"
#include "report/timer.hpp"

using namespace std;

PerformanceAnalyzer::PerformanceAnalyzer(boost::mpi::communicator thisworld, std::map<boost::thread::id,Timer>& timermap) {
	// aggregate all timer on head node
//	boost::mpi::gather(thisworld,timer,m_alltimer,0);

	// we need: sum , count, mean , variance, min, max
	double my_sum, my_count, my_mean, my_variance, my_min, my_max;
	vector<double> all_sum, all_count, all_mean, all_variance, all_min, all_max;

	int my_keyflag;
	vector<int> all_keyflags;
	
	//	vector<double> all_sum, all_count, all_mean, all_variance, all_min, all_max;
	all_sum.assign(thisworld.size(),0);
	all_count.assign(thisworld.size(),0);
	all_mean.assign(thisworld.size(),0);
	all_variance.assign(thisworld.size(),0);
	all_min.assign(thisworld.size(),0);
	all_max.assign(thisworld.size(),0);

	all_keyflags.assign(thisworld.size(),0);

	// first negotiate keys 
    set<string> my_keys_set;
    for(std::map<boost::thread::id,Timer>::iterator i = timermap.begin(); i != timermap.end(); ++i)
    {
        vector<string> keys = i->second.keys();
        for(size_t i = 0; i < keys.size(); ++i)
        {
            my_keys_set.insert(keys[i]);            
        }
    }
    vector<string> my_keys;
    for(set<string>::iterator i = my_keys_set.begin(); i != my_keys_set.end(); ++i)
    {
        my_keys.push_back(*i);
    }
    
	vector< vector<string> > all_keys(thisworld.size());

	boost::mpi::gather(thisworld ,my_keys,all_keys,0);
	
	set<string> common_keys_set;
	if (thisworld.rank()==0) {
		for(size_t i = 0; i < all_keys.size(); ++i)
		{	
			for(vector<string>::iterator j = all_keys[i].begin(); j != all_keys[i].end(); ++j)
			{
                common_keys_set.insert(*j);
			}
		}		
	}
	vector<string> common_keys;
    for(set<string>::iterator i = common_keys_set.begin(); i != common_keys_set.end(); ++i)
    {
        common_keys.push_back(*i);
    }
    mpi::wrapper::broadcast_class<vector<string> >(thisworld,common_keys,0);
		
	for(vector<string>::iterator ti = common_keys.begin(); ti != common_keys.end(); ++ti)
	{
		PerformanceMeasure& pm = m_measures[*ti];
		pm.sum = 0;
		pm.count = 0;
		pm.mean = 0;
		pm.variance = 0;
		pm.min = 0;
		pm.max = 0;

        // iterate over all timer per node first
        my_sum = 0; my_count = 0; my_mean = 0; my_variance = 0; my_min = 0; my_max = 0;
        my_keyflag = 0; 
    	double this_sum, this_count, this_mean, this_variance, this_min, this_max;
        size_t timer_count = 0;
        for(std::map<boost::thread::id,Timer>::iterator thistimer = timermap.begin(); thistimer != timermap.end(); ++thistimer)
        {
            Timer& timer = thistimer->second;
            vector<string> keys = timer.keys();
            size_t this_keyflag = 0;
            if (find(keys.begin(),keys.end(),*ti)!=keys.end()) {
                this_keyflag = 1;
            }

    		if (this_keyflag == 1) {
                timer_count++;
                
    			this_sum = timer.sum(*ti);
    			this_count = timer.count(*ti);
    			this_mean = timer.mean(*ti);
    			this_variance = timer.variance(*ti);
    			this_min = timer.min(*ti);
    			this_max = timer.max(*ti);			

                my_sum += this_sum;
                my_count += this_count;
                if (my_keyflag == 0 ) {
                    my_min = this_min;
                    my_max = this_max;
                } else {
                    if (my_min>this_min) my_min = this_min;
                    if (my_max<this_max) my_max = this_max;   
                }
                // this mean has to be weighted by the count
                // correction has to happen afterwards
                my_mean += this_mean*this_count;
                my_variance += this_variance*this_count;

                my_keyflag = 1;
    		}            
        }

		if (my_keyflag == 1) {
            my_mean /= my_count;
            my_variance /= my_count;
		}
		
		boost::mpi::gather(thisworld ,my_keyflag  , all_keyflags ,0);		
		boost::mpi::gather(thisworld ,my_sum      , all_sum      ,0);
		boost::mpi::gather(thisworld ,my_count    , all_count    ,0);
		boost::mpi::gather(thisworld ,my_mean     , all_mean     ,0);
		boost::mpi::gather(thisworld ,my_variance , all_variance ,0);
		boost::mpi::gather(thisworld ,my_min      , all_min      ,0);
		boost::mpi::gather(thisworld ,my_max      , all_max      ,0);

        size_t NN = thisworld.size();
		if (thisworld.rank()==0) {
			double weightedmean=0;
			double weightedvariance=0;
			double totalmin = 1e80;
			double totalmax = -1;

			for(size_t i = 0; i < NN; ++i)
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

// end of file
