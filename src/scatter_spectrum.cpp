/*
 *  scatterspectrum.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_spectrum.hpp"

// standard header
#include <fstream>

// special library headers

// other headers
#include "parameters.hpp"

using namespace std;

ScatterSpectrum::ScatterSpectrum(std::vector<ScatterSpectrum> scatspecs) {
	for(size_t i = 0; i < scatspecs.size(); ++i)
	{
		for(size_t j = 0; j < scatspecs[i].size(); ++j)
		{
			this->push_back(scatspecs[i][j]);
		}
	}
//	sort(this->begin(),this->end());
}


void ScatterSpectrum::write_plain(string fname,string format) {
	
	ofstream ofile(fname.c_str());

	double supermaxvalue = 0;
	std::map<CartesianCoor3D,double> maxvalues;
	// first scan for the maximum value, set it to 1.0
	for (ScatterSpectrum::iterator si=this->begin();si!=this->end();si++) {
		const CartesianCoor3D& qvector = si->first;
		double maxvalue = 0;
		for (std::vector<std::complex<double> >::iterator si2=si->second.begin();si2!=si->second.end();si2++) {
			const complex<double>& value = *si2;
			if (value.real()>maxvalue) maxvalue = value.real();		
		}
		maxvalues[qvector] = maxvalue;
		if (maxvalue>supermaxvalue) supermaxvalue = maxvalue;
	}
	
	
	if (format == "txt") {
		// normalized value: I(frame)/I_max_overall/I_max_qvector
		ofile << "#" << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "frame" << "\t" << "I_max_overall" << "\t" << "I_max_qvector" << "\t" << "I(frame)" << endl;
	}
	else {
		throw;
	}


	for (ScatterSpectrum::iterator si=this->begin();si!=this->end();si++) {
		// review this... declare a prooper "trajectory format!"
		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
		const CartesianCoor3D& qvector = si->first;

		size_t framenumber = 0;
		for (std::vector<std::complex<double> >::iterator si2=si->second.begin();si2!=si->second.end();si2++) {
			const double& value = si2->real();

			std::vector<std::complex<double> >::iterator testiterator = si2;			
			
			if (format == "txt") {
				ofile <<  qvector.x << "\t" <<  qvector.y << "\t" <<  qvector.z <<  "\t" << framenumber << "\t" << supermaxvalue << "\t" << maxvalues[qvector] << "\t"  << value <<  endl;			
				if ( testiterator++ == si->second.end()) { ofile << endl; } // if we hit the end of frames for current vector...
			}
			
			framenumber++;
		}
		ofile << endl;
	}
	
};



void ScatterSpectrum::write_average(string fname,string format) {
	
	ofstream ofile(fname.c_str());

	// if this run was made via the scan feature, include some gaps for easy processing in gnuplot
	size_t scansize1,scansize2,scansize3;
	scansize1 = scansize2 = scansize3 = this->size();
	if (Params::Inst()->scattering.qvectors.scans.size()>0) {
		std::vector<ScatteringVectorsScanParameters>& scans = Params::Inst()->scattering.qvectors.scans;
		if (scans.size()==1) {
			scansize1 = scans[0].points;
		}
		if (scans.size()==2) {
			scansize1 = scans[1].points;
			scansize2 = scans[0].points;
		}
		if (scans.size()==3) {
			scansize1 = scans[2].points;
			scansize2 = scans[1].points;
			scansize3 = scans[0].points;						
		}

	}

	double max = 0;	
	vector<pair<CartesianCoor3D,pair<double,double> > > tempvalues;
	// first scan for the maximum value, set it to 1.0
	for (ScatterSpectrum::iterator si=this->begin();si!=this->end();si++) {
		
		const CartesianCoor3D& qvector = si->first;
		vector<complex<double> >& values = si->second;
		double sum=0;
		double sum2=0;
		for(size_t i = 0; i < values.size(); ++i)
		{
			sum += values[i].real();
			sum2 += values[i].real()*values[i].real();
		}

		double mean = sum/values.size();
		double var = sum2/values.size() - mean*mean;
		if (mean>max) max = mean;
		tempvalues.push_back(make_pair(qvector,make_pair(mean,var)));
	}
	
	
	if (format == "txt") {
		// normalized value: I(frame)/I_max_overall/I_max_qvector
		ofile << "#" << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "I_max" << "\t" << "I(frame)" << "\t" << "I(Stddev)" << endl;
	}
	else {
		throw;
	}

	size_t count = 1;
	for(size_t i = 0; i < tempvalues.size(); ++i)
	{
		tempvalues[i].first;
		
		if (format == "txt") {
			ofile <<  tempvalues[i].first.x 			<< "\t";
			ofile <<  tempvalues[i].first.y 			<< "\t";
			ofile <<  tempvalues[i].first.z 			<< "\t";
			ofile <<  max 								<< "\t";
			ofile <<  tempvalues[i].second.first 		<< "\t";
			ofile <<  sqrt(tempvalues[i].second.second) << "\t";
			ofile << endl;			
		}
		if (count==tempvalues.size()) break;

		if ( (count % (scansize1))==0) ofile << endl;
		if ( (count % (scansize1*scansize2))==0) ofile << endl;
		if ( (count % (scansize1*scansize2*scansize3))==0) ofile << endl;

		count++;
	}
};

// end of file
