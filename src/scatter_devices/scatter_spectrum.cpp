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
#include "scatter_devices/scatter_spectrum.hpp"

// standard header
#include <fstream>

// special library headers
#include <fftw3.h>
#include <boost/algorithm/string.hpp>
#include <H5Cpp.h>

// other headers
#include "control.hpp"
#include "log.hpp"

using namespace std;

ScatterSpectrum::ScatterSpectrum(std::vector<ScatterSpectrum> scatspecs) {
	for(size_t i = 0; i < scatspecs.size(); ++i)
	{
		for(size_t j = 0; j < scatspecs[i].size(); ++j)
		{
            qvectors.push_back(scatspecs[i].qvectors[j]);
            amplitudes.push_back(scatspecs[i].amplitudes[j]);
		}
	}
//	sort(this->begin(),this->end());
}

void ScatterSpectrum::add(CartesianCoor3D q, std::vector<std::complex<double> > spectrum) {
    qvectors.push_back(q);
    amplitudes.push_back(spectrum);
}

void ScatterSpectrum::clear() {
    qvectors.clear();
    amplitudes.clear();
}

void ScatterSpectrum::read_plain(string fname,string format) {
   
   clear();
    
	if (format == "txt") {

        std::map<CartesianCoor3D,std::vector<std::complex<double> > > values;
	    
	   ifstream ifile(fname.c_str());
        
       string line;
	    // skip empty lines and lines starting with #
       while (getline(ifile,line)) {
           boost::trim(line);
           if (line=="") continue;
           if (boost::starts_with(line,"#")) continue;
           double qx,qy,qz,frame,Imax_overall,Imax_qvector,I0_qvector,Ir,Ii;
           stringstream linestream(line);
           linestream >> qx >> qy >> qz >> frame >> Imax_overall >> Imax_qvector >> I0_qvector >> Ir >> Ii;
           values[CartesianCoor3D(qx,qy,qz)].push_back(std::complex<double>(Ir,Ii));
       }
       
       for(std::map<CartesianCoor3D,std::vector<std::complex<double> > >::iterator vi=values.begin();vi!=values.end();vi++) {
           add(vi->first,vi->second);
       }
   } else if (format == "hdf5") {
       
   }
}

void ScatterSpectrum::write_plain(string fname,string format) {

	double supermaxvalue = 0;
	std::map<CartesianCoor3D,double> maxvalues;
	// first scan for the maximum value, set it to 1.0
	for(size_t di = 0; di < qvectors.size(); ++di)
	{
		const CartesianCoor3D& qvector = qvectors[di];
		double maxvalue = 0;
		for (std::vector<std::complex<double> >::iterator si2=amplitudes[di].begin();si2!=amplitudes[di].end();si2++) {
			const complex<double>& value = *si2;
			if (value.real()>maxvalue) maxvalue = value.real();		
		}
		maxvalues[qvector] = maxvalue;
		if (maxvalue>supermaxvalue) supermaxvalue = maxvalue;
	}
	
	if (format == "txt") {
	    
	    ofstream ofile(fname.c_str());
    	
		// normalized value: I(frame)/I_max_overall/I_max_qvector
		ofile << "#" << "qx" << "\t" << "qy" << "\t" << "qz" << "\t" << "frame" << "\t" << "I_max_overall" << "\t" << "I_max_qvector" <<  "\t" << "I_0_qvector" << "\t" << "I(frame)" << endl;
		
    	for(size_t di = 0; di < qvectors.size(); ++di)
    	{
    		// review this... declare a prooper "trajectory format!"
    		// espeically clarify whether complex number should be written and how the standard deviations are interpreted
    		const CartesianCoor3D& qvector = qvectors[di];

    		size_t framenumber = 0;
    		complex<double> Iq0;
    		if (size()>0) Iq0 = amplitudes[di].begin()->real();
    		for (std::vector<std::complex<double> >::iterator si2=amplitudes[di].begin();si2!=amplitudes[di].end();si2++) {

    			std::vector<std::complex<double> >::iterator testiterator = si2;			

    			if (format == "txt") {
    				ofile <<  qvector.x << "\t" <<  qvector.y << "\t" <<  qvector.z <<  "\t" << framenumber << "\t" << supermaxvalue << "\t" << maxvalues[qvector] << "\t" << Iq0.real() << "\t"  << si2->real() << "\t"  << si2->imag()  <<  endl;			
    				if ( testiterator++ == amplitudes[di].end()) { ofile << endl; } // if we hit the end of frames for current vector...
    			}

    			framenumber++;
    		}
    		ofile << endl;
    	}		
	}
    
	else if (format == "hdf5") {

        // create header
        H5::H5File h5file( fname, H5F_ACC_TRUNC );

        hsize_t qvector_field[2];
        qvector_field[0]=qvectors.size();
        qvector_field[1]=3;
        H5::DataSpace dspace_qvectors(2,qvector_field);

        
        H5::DataType dtype_qvectors(H5::PredType::NATIVE_DOUBLE);
	    H5::DataSet dset_qvectors = h5file.createDataSet( "qvectors", dtype_qvectors, dspace_qvectors );
        dset_qvectors.write(reinterpret_cast<double*>(&qvectors[0]), H5::PredType::NATIVE_DOUBLE);


        // detect the number of frames from first amplitudes
        size_t NF = 0;
        if (qvectors.size()>0) NF = amplitudes[0].size();
        
        hsize_t fqt_field[3];
        fqt_field[0]=qvectors.size();
        fqt_field[1]=NF;
        fqt_field[2]=2;
        H5::DataSpace dspace_fqt(3,fqt_field);
        
        H5::DataType dtype_fqt(H5::PredType::NATIVE_DOUBLE);
	    H5::DataSet dset_fqt = h5file.createDataSet( "fqt", dtype_fqt, dspace_fqt );

        hsize_t start[3];  // Start of hyperslab
        hsize_t stride[3]; // Stride of hyperslab
        hsize_t count[3];  // Block count
        hsize_t block[3];  // Block sizes
        start[0]=0;start[1]=0;start[2]=0;
        count[0]=1;count[1]=NF;count[2]=2;
        
        hsize_t mfqt_field[1];
        mfqt_field[0] = 2*NF;
        H5::DataSpace mspace_fqt(1,mfqt_field);

        for(size_t i=0;i<qvectors.size();i++) {
            start[0]=i;
            dspace_fqt.selectNone();
            dspace_fqt.selectHyperslab(H5S_SELECT_SET,count,start);
            
            dset_fqt.write(reinterpret_cast<double*>(&amplitudes[i][0].real()), H5::PredType::NATIVE_DOUBLE,mspace_fqt, dspace_fqt);                        
        }
        
	}
	else {
		throw;
	}




	
};


void ScatterSpectrum::write_average(string fname,string format) {
	
	ofstream ofile(fname.c_str());

	// if this run was made via the scan feature, include some gaps for easy processing in gnuplot
	size_t scansize1,scansize2,scansize3;
	scansize1 = scansize2 = scansize3 = size();
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
	for(size_t di = 0; di < size(); ++di)
	{
		
		const CartesianCoor3D& qvector = qvectors[di];
		vector<complex<double> >& values = amplitudes[di];
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
