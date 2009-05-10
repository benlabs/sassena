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

//bool ScatterdataKey::operator<(const ScatterdataKey that) const{
//	// put lower frame first
//	if (this->second < that.second)
//		return true;
//
//	if (this->second > that.second)
//		return false;
//
//	// then cluster by vectors
//	// x being first iso line
//	if (this->first.x < that.first.x)
//		return true;
//
//	if (this->first.x > that.first.x)
//		return false;
//
//	if (this->first.y < that.first.y)
//		return true;
//
//	if (this->first.y > that.first.y)
//		return false;
//
//	if (this->first.z < that.first.z)
//		return true;
//
//	if (this->first.z > that.first.z)
//		return false;
//
//
//	return false;
//}

template <class V> void Scatterdata<V>::write_plain(string fname,string format) {
	
	ofstream ofile(fname.c_str());

	V supermaxvalue = 0;
	std::map<CartesianCoor3D,V> maxvalues;
	// first scan for the maximum value, set it to 1.0
	for (typename Scatterdata<V>::iterator si=this->begin();si!=this->end();si++) {
		const CartesianCoor3D& qvector = si->first;
		V maxvalue = 0;
		for (typename std::map<size_t,V>::iterator si2=si->second.begin();si2!=si->second.end();si2++) {
			const V& value = si2->second;
			size_t frame = si2->first;
			if (value>maxvalue) maxvalue = value;		
		}
		maxvalues[qvector] = maxvalue;
		if (maxvalue>supermaxvalue) supermaxvalue = maxvalue;
	}
	
	
	if (format == "txt") {
		ofile << "#" << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "frame" << "\t" << "I_max_overall" << "\t" << "I_max_qvector" << "\t" << "I(frame)/I_max_overall/I_max_qvector" << endl;
	}
	else {
		throw;
	}


	for (typename Scatterdata<V>::iterator si=this->begin();si!=this->end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const CartesianCoor3D& qvector = si->first;
		
		for (typename std::map<size_t,V>::iterator si2=si->second.begin();si2!=si->second.end();si2++) {
			const V& value = si2->second;
			size_t frame = si2->first;
			typename std::map<size_t,V>::iterator testiterator = si2;			
			
			if (format == "txt") {
				ofile <<  qvector.x << "\t" <<  qvector.y << "\t" <<  qvector.z <<  "\t" << frame << "\t" << supermaxvalue << "\t" << maxvalues[qvector] << "\t"  << value/maxvalues[qvector]/supermaxvalue <<  endl;			
				if ( testiterator++ == si->second.end()) { ofile << endl; } // if we hit the end of frames for current vector...
			}
		}
		ofile << endl;
	}
	
};


template <class V> void Scatterdata<V>::write_average(string fname,string format) {
	
	ofstream ofile(fname.c_str());
	
	// use helper method average
	typedef std::map<CartesianCoor3D,std::pair<V,V> > mappair;
	mappair	a = average();

	V maxvalue=0;
	// first scan for the maximum value, set it to 1.0
	for (typename mappair::iterator si=a.begin();si!=a.end();si++) {
		const pair<V,V>& value = si->second;
		if (value.first>maxvalue) maxvalue = value.first;		
	}

	if (format == "txt") {
		ofile << "#" << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "I_max" << "\t" << "I/I_max" << "\t" << "stddev(I)/I_max" << endl;
	}
	else {
		throw;
	}

	for (typename mappair::iterator si=a.begin();si!=a.end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const CartesianCoor3D& qvector = si->first;
		const pair<V,V>& value = si->second;
		
		if (format == "txt") {
			ofile <<  qvector.x << "\t" <<  qvector.y << "\t" <<  qvector.z << "\t" << maxvalue <<  "\t" << value.first/maxvalue << "\t" << value.second/maxvalue <<  endl;			
		}
	}
	
};

template <class V> std::map<CartesianCoor3D,std::pair<V,V> > Scatterdata<V>::average() {

	using namespace boost::accumulators;
			
	typedef std::map<CartesianCoor3D,accumulator_set<V,features<tag::mean,tag::variance> > > mapacc;
	mapacc r;

	for (typename Scatterdata<V>::iterator sdi=this->begin();sdi!=this->end();sdi++) {
		for (typename map<size_t,V>::iterator sdi2=sdi->second.begin();sdi2!=sdi->second.end();sdi2++) {		
			r[sdi->first](sdi2->second); // calls the operator() of the accumulator
		}
	}

	// return type:
	std::map<CartesianCoor3D,pair<V,V> > ret;
	for (typename mapacc::iterator ri=r.begin();ri!=r.end();ri++) {
		ret[ri->first]=make_pair(mean(ri->second),sqrt(variance(ri->second)));
	}
	
	return ret;
}

template <class V> void Scatterdata<V>::dump(streambuf* pof) {

	ostream os(pof);
	
	os << "INFO>> " << "FORMAT: " << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "frame" << "\t" << "value" << endl;

	for (typename Scatterdata<V>::iterator si=this->begin();si!=this->end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const CartesianCoor3D& qvector = si->first;
		for (typename std::map<size_t,V>::iterator si2=si->second.begin();si2!=si->second.end();si2++) {
			const V& value = si2->second;
			size_t frame = si2->first;
			os <<  qvector.x << "\t" <<  qvector.y << "\t" <<  qvector.z << "\t" << frame << "\t" << value <<  endl;
		}
	}
	
};

// explicit instantiation
template class Scatterdata<double>;
//template class Scatterdata<complex<double> >;

// end of file


