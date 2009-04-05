/*
 *  scatterdata.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatterdata.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers

using namespace std;

bool ScatterdataKey::operator<(const ScatterdataKey that) const{
	// put lower frame first
	if (this->second < that.second)
		return true;

	if (this->second > that.second)
		return false;

	// then cluster by vectors
	// x being first iso line
	if (this->first.x < that.first.x)
		return true;

	if (this->first.x > that.first.x)
		return false;

	if (this->first.y < that.first.y)
		return true;

	if (this->first.y > that.first.y)
		return false;

	if (this->first.z < that.first.z)
		return true;

	if (this->first.z > that.first.z)
		return false;


	return false;
}

template <class V> void Scatterdata<V>::write_plain(string fname,string format) {
	
	ofstream ofile(fname.c_str());
	
	for (typename Scatterdata<V>::iterator si=this->begin();si!=this->end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const ScatterdataKey& key = si->first;
		const V& value = si->second;
		
		if (format == "txt") {
			ofile <<  key.first.x << "\t" <<  key.first.y << "\t" <<  key.first.z << "\t" << key.second << "\t" << value <<  endl;			
		}
		else {
			cerr << "ERROR>> " << "output format not understood" << endl;
		}
	}
	
};


template <class V> void Scatterdata<V>::write_average(string fname,string format) {
	
	ofstream ofile(fname.c_str());
	
	// use helper method average
	typedef std::map<ScatterdataKey,std::pair<V,V> > mappair;
	mappair	a = average();

	for (typename mappair::iterator si=a.begin();si!=a.end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const ScatterdataKey& key = si->first;
		const pair<V,V>& value = si->second;
		
		if (format == "txt") {
			ofile <<  key.first.x << "\t" <<  key.first.y << "\t" <<  key.first.z << "\t" << value.first << "\t" << value.second <<  endl;			
		}
		else {
			cerr << "ERROR>> " << "output format not understood" << endl;
		}
	}
	
};

template <class V> std::map<ScatterdataKey,std::pair<V,V> > Scatterdata<V>::average() {

	using namespace boost::accumulators;
			
	typedef std::map<ScatterdataKey,accumulator_set<V,features<tag::mean,tag::variance> > > mapacc;
	mapacc r;

	for (typename Scatterdata<V>::iterator sdi=this->begin();sdi!=this->end();sdi++) {
		r[ScatterdataKey(sdi->first.first,0)](sdi->second); // calls the operator() of [key(coor,0)] accumulator
	}

	// return type:
	std::map<ScatterdataKey,pair<V,V> > ret;
	for (typename mapacc::iterator ri=r.begin();ri!=r.end();ri++) {
		ret[ri->first]=make_pair(mean(ri->second),sqrt(variance(ri->second)));
	}
	
	return ret;
}


template <class V> void Scatterdata<V>::dump(streambuf* pof) {

	ostream os(pof);

	for (typename Scatterdata<V>::iterator si=this->begin();si!=this->end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const ScatterdataKey& key = si->first;
		const V& value = si->second;
		os <<  key.first.x << "\t" <<  key.first.y << "\t" <<  key.first.z << "\t" << key.second << "\t" << value <<  endl;
	}
	
};

// explicit instantiation
template class Scatterdata<double>;
//template class Scatterdata<complex<double> >;

// end of file


