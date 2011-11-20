/** \file
This file contains an efficient timer class, which is used to retrieve execution times for various parts of the algorithms.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef REPORT__TIMER_HPP_
#define REPORT__TIMER_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers

using namespace boost::accumulators;

/** 
Type class which respresents the time value used by the Timer class
*/
class Timer_timeval {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & tv_sec;
		ar & tv_usec;
    }
	/////////////////// 
public:
	long tv_sec;
	long tv_usec;
	
	Timer_timeval() {}
	Timer_timeval(const timeval& t) : tv_sec(t.tv_sec) , tv_usec(t.tv_usec) {}
	Timer_timeval& operator=(const timeval& t) { tv_sec = t.tv_sec; tv_usec = t.tv_usec; return (*this); }
};

/** 
Basic Timer which provides a start/stop facility measure runtimes. Provides an interface to retrieve statistical information.
*/
class Timer {
	private:
		/////////////////// MPI related
		// make this class serializable to 
		// allow sample to be transmitted via MPI
	    friend class boost::serialization::access;	
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
	    {
			ar & times;
			ar & starttimes;
			ar & states;
	    }
		/////////////////// 
	
	typedef boost::accumulators::accumulator_set<double,features<tag::min,tag::max,tag::mean,tag::variance,tag::sum,tag::count> > times_type;
	std::map<std::string,times_type > times;
//	std::map<std::string,std::vector<double> > times;
	
	std::map<std::string,Timer_timeval> starttimes;
	std::map<std::string,bool> states;
	
	double t_diff(timeval start, timeval end);
	double t_diff(Timer_timeval start, timeval end);

public:			
	std::vector<std::string> keys();
	
	void start(std::string tk);
	void stop(std::string tk);

	double sum(std::string tk);
	double min(std::string tk);
	double max(std::string tk);	
	double mean(std::string tk);
	double variance(std::string tk);	
	double count(std::string tk);
	
	void clear();
	
	bool has_key(std::string tk);
//	Timer& operator+=( Timer& other);
//	Timer operator+( Timer& other);
};

#endif 

// end of file
