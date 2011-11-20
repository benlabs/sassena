/** \file
This file contains mathemical routines.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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


void auto_correlate_direct(fftw_complex* data,size_t N) {

    size_t NF = N;

    fftw_complex* data_local = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
    memcpy(data_local,data,N*sizeof(fftw_complex));

    // direct
    for(size_t tau = 0; tau < NF; ++tau)
    {
        data[tau][0]=0;data[tau][1]=0;
    	size_t last_starting_frame = NF-tau;
    	for(size_t k = 0; k < last_starting_frame; ++k)
    	{
    		fftw_complex& a1 = data_local[k];
    		fftw_complex& a2 = data_local[k+tau];
            data[tau][0] += a1[0]*a2[0]+a1[1]*a2[1];
            data[tau][1] += -a1[0]*a2[1]+a1[1]*a2[0];            
    	}
    	data[tau][0] /= (last_starting_frame); 		
    	data[tau][1] /= (last_starting_frame); 		    	
    }    
    
    fftw_free(data_local);

}

void auto_correlate_fftw(std::vector<std::complex<double> >& data,size_t N,fftw_plan p1,fftw_plan p2,fftw_complex* fftw_planspace) {
    size_t NF = data.size();

    memcpy(fftw_planspace,&(data[0]),sizeof(double)*2*NF);
    for(size_t i = NF; i < 2*NF; ++i)
    {
        fftw_planspace[i][0]=0;
        fftw_planspace[i][1]=0;        
    }

    fftw_execute_dft(p1,fftw_planspace,fftw_planspace);
    for(size_t i = 0; i < 2*NF; ++i)  {
        fftw_planspace[i][0]=fftw_planspace[i][0]*fftw_planspace[i][0]+fftw_planspace[i][1]*fftw_planspace[i][1];
        fftw_planspace[i][1]=0;  
    }
    fftw_execute_dft(p2,fftw_planspace,fftw_planspace);

    for(size_t i = 0; i < NF; ++i) {
        double norm = ( 2*NF * (NF -i ) );
        data[i] = std::complex<double>(fftw_planspace[i][0]/norm,fftw_planspace[i][1]/norm); 
    }    
}


template <class T> void auto_correlate_fftw(std::vector<std::complex<T> >& data,fftw_plan p1,fftw_plan p2) {
    size_t NF = data.size();
    
    // make a local copy to allow to override data
    std::vector<std::complex<T> > data_local = data;
    
    std::vector< std::complex<T> >& complete_a = data_local;
    std::vector< std::complex<T> >& correlated_a = data;
    
    
    fftw_complex *wspace;
    wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*NF);


    for(size_t i = 0; i < NF; ++i) {
        wspace[i][0]= complete_a[i].real();            
        wspace[i][1]= complete_a[i].imag();                        
    }
    for(size_t i = NF; i < 2*NF; ++i) {
        wspace[i][0]=0;
        wspace[i][1]=0;
    }
    
    fftw_execute_dft(p1,wspace,wspace); /* repeat as needed */
    for(size_t i = 0; i < 2*NF; ++i)  {
        wspace[i][0]=wspace[i][0]*wspace[i][0]+wspace[i][1]*wspace[i][1];
        wspace[i][1]=0;  
    }
    fftw_execute_dft(p2,wspace,wspace); /* repeat as needed */

    for(size_t i = 0; i < NF; ++i) {
        correlated_a[i]=std::complex<T>(wspace[i][0],wspace[i][1])*( T(1.0) / ( 2*NF * (NF -i ) ) );
    }
    
    fftw_free(wspace);
    
}


void auto_correlate_fftw(fftw_complex* data,fftw_plan p1,fftw_plan p2,size_t NF) {
        
    fftw_execute_dft(p1,data,data); /* repeat as needed */
    for(size_t i = 0; i < 2*NF; ++i)  {
        data[i][0]=data[i][0]*data[i][0]+data[i][1]*data[i][1];
        data[i][1]=0;  
    }
    fftw_execute_dft(p2,data,data); /* repeat as needed */

    for(size_t i = 0; i < NF; ++i) {
        double factor = (1.0/(2*NF*(NF-i)));
        data[i][0]*=factor;
        data[i][1]*=factor;        
    }
    
}

template <class T> void square_elements(std::vector<std::complex<T> >& data) {

    size_t NF = data.size();
    for(size_t n = 0; n < NF; n++)
    {
        data[n]*=conj(data[n]);
    }
}


void square_elements(fftw_complex* data,size_t N) {

    size_t NF = N;
    for(size_t n = 0; n < NF; n++)
    {
        double r = data[n][0]*data[n][0]+data[n][1]*data[n][1];
        data[n][0]=r;
        data[n][1]=0;        
    }
}


template <class T> void multiply_elements(const T& factor,std::vector<std::complex<T> >& data) {

    size_t NF = data.size();
    for(size_t n = 0; n < NF; n++)
    {
        data[n]*=factor;
    }
}


void multiply_elements(const double factor,fftw_complex* data,size_t NF) {
    for(size_t n = 0; n < NF; n++)
    {
        data[n][0]*=factor;
        data[n][1]*=factor;
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


void add_elements(fftw_complex* target,const fftw_complex* source,size_t NF) {
    for(size_t n = 0; n < NF; n++)
    {
        target[n][0]+=source[n][0];
        target[n][1]+=source[n][1];
    }
}

template<class T> void add_elements(std::vector<std::complex<T> >& target,const fftw_complex* source,size_t N) {
    
    size_t NF = target.size();
    if (N<NF) NF=N;
    for(size_t n = 0; n < NF; n++)
    {
        target[n]+=std::complex<T>(source[n][0],source[n][1]);
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


template <class T> std::complex<T> reduce(const fftw_complex* data,size_t N) {
    size_t size = N;
    std::complex<T> s=0;
    for(size_t i = 0; i < size; ++i)
    {
        s+=std::complex<T>(data[i][0],data[i][1]);
    }
    return s;
}

}



template void smath::square_elements<float>(std::vector<std::complex<float> >& data);
template void smath::square_elements<double>(std::vector<std::complex<double> >& data);

template void smath::multiply_elements<float>(const float& factor,std::vector<std::complex<float> >& data);
template void smath::multiply_elements<double>(const double& factor,std::vector<std::complex<double> >& data);

template void smath::add_elements<float>(std::vector<std::complex<float> >& target,const std::vector<std::complex<float> >& source);
template void smath::add_elements<double>(std::vector<std::complex<double> >& target,const std::vector<std::complex<double> >& source);

template void smath::add_elements<float>(std::vector<std::complex<float> >& target,const fftw_complex* source,size_t N);
template void smath::add_elements<double>(std::vector<std::complex<double> >& target,const fftw_complex* source,size_t N);

template void smath::auto_correlate_fftw<float>(std::vector<std::complex<float> >& data,fftw_plan p1,fftw_plan p2);
template void smath::auto_correlate_fftw<double>(std::vector<std::complex<double> >& data,fftw_plan p1,fftw_plan p2);
              
template void smath::auto_correlate_direct<float>(std::vector<std::complex<float> >& data);
template void smath::auto_correlate_direct<double>(std::vector<std::complex<double> >& data);
              
template std::complex<float> smath::reduce<float>(const std::vector<std::complex<float> >& data);
template std::complex<double> smath::reduce<double>(const std::vector<std::complex<double> >& data);

template std::complex<float> smath::reduce<float>(const fftw_complex* data,size_t N);
template std::complex<double> smath::reduce<double>(const fftw_complex* data,size_t N);


// end of file
