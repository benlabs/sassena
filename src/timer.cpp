#include "timer.hpp"

#include <deque>

#include "log.hpp"
#include "parameters.hpp"

using namespace std;

double Timer::t_diff(timeval start, timeval end) {
	int seconds = (end.tv_sec-start.tv_sec);
	int microseconds = (end.tv_usec-start.tv_usec);
	seconds = (microseconds < 0) ? seconds-1 : seconds;
	microseconds = (microseconds < 0) ? 1000000+microseconds : microseconds;
	return seconds + (microseconds/1000000.0);
}

void Timer::start(std::string tk)  { 
	states[tk]=true;

	timeval t;
	gettimeofday(&t, 0);
	starttimes[tk] = t;
	
	if (Params::Inst()->debug.timer) Info::Inst()->write(string("Starting timer for <")+tk+string(">"));
}

void Timer::stop(std::string tk)   { 
	if (states.find(tk)==states.end()) return; 
	if (states[tk]) { 
		states[tk]=false;
		timeval t;
		gettimeofday(&t, 0);
		times[tk](t_diff(starttimes[tk],t));
	}
	else {
		cerr << "ERROR>> " << " Timer not started, but stopped: " << tk << endl;
		throw;
	}
	if (Params::Inst()->debug.timer) Info::Inst()->write(string("Stopping timer for <")+tk+string(">"));
}

vector<string> Timer::keys() {
	vector<string> ret;
	for (map<string,times_type>::iterator ti=times.begin();ti!=times.end();ti++) {
		ret.push_back(ti->first);
	}
	return ret;
}

double Timer::mean(std::string tk) {
	return boost::accumulators::mean(times[tk]);
}

double Timer::variance(std::string tk) {
	return boost::accumulators::variance(times[tk]);
}

double Timer::sum(std::string tk) {
	return boost::accumulators::sum(times[tk]);
}
double Timer::min(std::string tk) {
	return boost::accumulators::min(times[tk]);
}
double Timer::max(std::string tk) {
	return boost::accumulators::max(times[tk]);
}

double Timer::count(std::string tk) {
	return boost::accumulators::count(times[tk]);
}

// end of file
