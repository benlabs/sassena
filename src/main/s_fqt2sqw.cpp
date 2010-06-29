/*
 *  sassena.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "common.hpp"

// standard header
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/math/distributions/normal.hpp>
#include <fftw3.h>
// other headers
#include "log.hpp"
#include "scatter_devices/scatter_spectrum.hpp"

using namespace std;

int main(int argc,char** argv) {

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(string("Info>>"));
	Warn::Inst()->set_prefix(string("Warn>>"));
	 Err::Inst()->set_prefix(string("Err>>"));
	
	if (argc!=4) {
        Err::Inst()->write("Need 2 arguments: FQT data file, energy resolution, boolean flag");
        throw;
	}

    ScatterSpectrum ss;
    
    
    ss.read_plain(argv[1],"txt");
    
    if (ss.size()<1) return 0;
    
    size_t seqlen = ss[0].second.size();
    
    // check sequence length
    for(size_t si = 0; si < ss.size(); ++si)
    {
        if (seqlen!=ss[si].second.size()) {
            Err::Inst()->write("Different Time Sequence Lengths. Exiting");
            return 1;
        }
    }
    
    size_t NF = seqlen;
    
 //   double dE = boost::lexical_cast<double>(argv[2]);
     double sigma_t = boost::lexical_cast<double>(argv[2]);
//    double sigma_e = 4*dE/(3*M_PI);
//    double sigma_t = 2/(M_PI*sigma_e);
    boost::math::normal normdist(0,sigma_t);
    vector<double> rt;

    // expand rt on demand
    if (NF>rt.size()) {
        size_t rtsize = rt.size();
        for(size_t i = rtsize; i < NF-rtsize; ++i) { 
            rt.push_back(boost::math::pdf(normdist,i)*sqrt(2*M_PI)*sigma_t);
        }
    }


    bool periodic = boost::lexical_cast<bool>(argv[3]);

    if (periodic) {

        fftw_complex *wspace;
        fftw_plan p1;
        wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NF-2));
        p1 = fftw_plan_dft_1d(2*NF-2, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    
        for(size_t si = 0; si < ss.size(); ++si)
        {
            Info::Inst()->write(string("qvector ")+to_s(si));
            vector<complex<double> >& fqt = ss[si].second;
                                    
            for(size_t i = 0; i < NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[i].real()*rt[i]; 
                wspace[i][1]= fqt[i].imag()*rt[i];                            
            }

            for(size_t i = NF; i < 2*NF-2; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[2*NF-2-i].real()*rt[2*NF-2-i]; 
                wspace[i][1]= fqt[2*NF-2-i].imag()*rt[2*NF-2-i];                            
            }

            fftw_execute(p1); /* repeat as needed */
            for(size_t i = 0; i < NF; ++i) {
                fqt[i]=0.5*complex<double>(wspace[i][0],wspace[i][1])*( 1.0 / ( 2*NF-2 ) );
            }
        }
        
        fftw_destroy_plan(p1);
        fftw_free(wspace);

    } else {
        fftw_complex *wspace;
        fftw_plan p1;
        wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NF));
        p1 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for(size_t si = 0; si < ss.size(); ++si)
        {
            Info::Inst()->write(string("qvector ")+to_s(si));
            vector<complex<double> >& fqt = ss[si].second;
            
            for(size_t i = 0; i < NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[i].real()*rt[i]; 
                wspace[i][1]= fqt[i].imag()*rt[i];                            

            }

            for(size_t i = NF; i < 2*NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= 0; 
                wspace[i][1]= 0;                            

            }
            fftw_execute(p1); /* repeat as needed */
            for(size_t i = 0; i < NF; ++i) {
                fqt[i]=0.5*complex<double>(wspace[i][0],wspace[i][1])*( 1.0 / ( 2*NF ) );
            }

        }        
        
        fftw_destroy_plan(p1);
        fftw_free(wspace);

    }



    
    for(size_t si = 0; si < ss.size(); ++si) {
        vector<complex<double> >& fqt = ss[si].second;
        for(size_t i = 0; i < fqt.size(); ++i)
        {
            cout << fqt[i].real() << endl;
        }
    }

	return 0;
}

// end of file