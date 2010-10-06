/*
 *  smath.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "math/smath.hpp"

// standard header
#include <cmath>
#include <iostream>

// special library headers

// other headers

using namespace std;

namespace smath {
    
template <class T> void auto_correlate_direct(std::vector<std::complex<T> >& data) {
    
    size_t NF = data.size();
    
    // make a local copy to allow to override data
    std::vector<std::complex<T> > data_local = data;
    
    std::vector< std::complex<T> >& complete_a = data_local;
    std::vector< std::complex<T> >& correlated_a = data;
            
    // direct
    for(size_t tau = 0; tau < NF; ++tau)
    {
        correlated_a[tau] = 0;
    	size_t last_starting_frame = NF-tau;
    	for(size_t k = 0; k < last_starting_frame; ++k)
    	{
    		std::complex<T>& a1 = complete_a[k];
    		std::complex<T>& a2 = complete_a[k+tau];
    		correlated_a[tau] += conj(a1)*a2;
    	}
    	correlated_a[tau] /= (last_starting_frame); 		
    }    
        
}

template <class T> void auto_correlate_fftw(std::vector<std::complex<T> >& data) {
    size_t NF = data.size();
    
    // make a local copy to allow to override data
    std::vector<std::complex<T> > data_local = data;
    
    std::vector< std::complex<T> >& complete_a = data_local;
    std::vector< std::complex<T> >& correlated_a = data;
    
    
    fftw_complex *wspace;
    fftw_plan p1,p2;
    wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*NF);
    p1 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);


    for(size_t i = 0; i < NF; ++i) {
        wspace[i][0]= complete_a[i].real();            
        wspace[i][1]= complete_a[i].imag();                        
    }
    for(size_t i = NF; i < 2*NF; ++i) {
        wspace[i][0]=0;
        wspace[i][1]=0;
    }
    
    fftw_execute(p1); /* repeat as needed */
    for(size_t i = 0; i < 2*NF; ++i)  {
        wspace[i][0]=wspace[i][0]*wspace[i][0]+wspace[i][1]*wspace[i][1];
        wspace[i][1]=0;  
    }
    fftw_execute(p2); /* repeat as needed */

    for(size_t i = 0; i < NF; ++i) {
        correlated_a[i]=std::complex<T>(wspace[i][0],wspace[i][1])*( T(1.0) / ( 2*NF * (NF -i ) ) );
    }
    
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(wspace);
    
}

template <class T> void square_elements(std::vector<std::complex<T> >& data) {

    size_t NF = data.size();
    for(size_t n = 0; n < NF; n++)
    {
        data[n]*=conj(data[n]);
    }
}


template <class T> void add_elements(std::vector<std::complex<T> >& target,const std::vector<std::complex<T> >& source) {
    
    size_t NF = target.size();
    if (source.size()<NF) NF=source.size();
    for(size_t n = 0; n < NF; n++)
    {
        target[n]+=source[n];
    }
}

template <class T> std::complex<T> reduce(const std::vector<std::complex<T> >& data) {
    size_t size = data.size();
    std::complex<T> s=0;
    for(size_t i = 0; i < size; ++i)
    {
        s+=data[i];
    }
    return s;
}

}



template void smath::square_elements<float>(std::vector<std::complex<float> >& data);
template void smath::square_elements<double>(std::vector<std::complex<double> >& data);

template void smath::add_elements<float>(std::vector<std::complex<float> >& target,const std::vector<std::complex<float> >& source);
template void smath::add_elements<double>(std::vector<std::complex<double> >& target,const std::vector<std::complex<double> >& source);

template void smath::auto_correlate_fftw<float>(std::vector<std::complex<float> >& data);
template void smath::auto_correlate_fftw<double>(std::vector<std::complex<double> >& data);
              
template void smath::auto_correlate_direct<float>(std::vector<std::complex<float> >& data);
template void smath::auto_correlate_direct<double>(std::vector<std::complex<double> >& data);
              
template std::complex<float> smath::reduce<float>(const std::vector<std::complex<float> >& data);
template std::complex<double> smath::reduce<double>(const std::vector<std::complex<double> >& data);
