/*
 *  scatterdevices.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "fftw/fftw++.h"
#include <fftw3.h>
#include "control.hpp"
#include "sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"
#include "vector_unfold.hpp"

using namespace std;

AllScatterDevice::AllScatterDevice(boost::mpi::communicator& thisworld, Sample& sample) {

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t nn = thisworld.size(); // Number of Nodes
	size_t na = sample.atoms.selections[target].size(); // Number of Atoms
	size_t nf = sample.coordinate_sets.size();

	size_t rank = thisworld.rank();

	EvenDecompose edecomp(nf,nn);
	
	myframes = edecomp.indexes_for(rank);

	if (thisworld.rank()==0) {	

        size_t memusage_scatmat = 2*sizeof(double)*myframes.size()*1;
                
        size_t memusage_per_cs = 3*sizeof(double)*na;
        size_t memusage_allcs = memusage_per_cs*myframes.size();
        
        
        Info::Inst()->write(string("Memory Requirements per node: "));
		Info::Inst()->write(string("Scattering Matrix: ")+to_s(memusage_scatmat)+string(" bytes"));


        // fault if not enough memory for scattering matrix
        if (memusage_scatmat>Params::Inst()->limits.memory.scattering_matrix) {
			Err::Inst()->write(string("Insufficient Buffer size for scattering matrix."));            
			Err::Inst()->write(string("Size required:")+to_s(memusage_scatmat)+string(" bytes"));            
			Err::Inst()->write(string("Configuration Parameter: limits.memory.scattering_matrix"));
            throw;
        }

		Info::Inst()->write(string("Coordinate Sets: ")+to_s(myframes.size()*memusage_per_cs)+string(" bytes"));		

        // warn if not enough memory for coordinate sets (cacheable system)
		if (Params::Inst()->runtime.limits.cache.coordinate_sets<myframes.size()) {
			Warn::Inst()->write(string("Insufficient Buffer size for coordinate sets. This is a HUGE bottleneck for performance!"));
			Warn::Inst()->write(string("Need at least: ")+to_s(memusage_allcs)+string(" bytes"));
			Warn::Inst()->write(string("Configuration Parameter: limits.memory.coordinate_sets"));
		}		
	}
	a.resize(myframes.size(),1); // n x 1 = frames x 1

	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
}

// acts like scatter_frame, but does the summation in place
void AllScatterDevice::scatter_frame_norm1(size_t localframe, CartesianCoor3D& q) {

	size_t noa = p_sample->coordinate_sets.get_selection().size();
	
	timer.start("sd:fs:f:ld");												
	CoordinateSet& cs = p_sample->coordinate_sets.load(myframes[localframe]);
	timer.stop("sd:fs:f:ld");
	
	vector<double>& sfs = scatterfactors.get_all();

	//	for(size_t j = 0; j < a.size2(); ++j)
	//	{
	//		double x1 = cs.x[j];
	//		double y1 = cs.y[j];
	//		double z1 = cs.z[j];	
	//		double p1 = x1*q.x+y1*q.y+z1*q.z;
	//
	//		a(iframe,j) = exp(-1.0*complex<double>(0,p1))*sfs[j];
	//	}



		double Ar = 0.0;
		double Ai = 0.0;
		double p, cp, ap;
		double M_PI_half = M_PI/2;
		double M_PI_3half = 3*M_PI/2;	
		double M_PI_twice = 2*M_PI;

		double csx = cs.c1[0];
		double csy = cs.c2[0];
		double csz = cs.c3[0];


		for(size_t j = 0; j < noa; ++j) {
			p =  csx*q.x + csy*q.y + csz*q.z;
//			double sign_sin = (p<0) ? -1.0 : 1.0;
//			ap = abs(p);
//			p = ap - long(ap/M_PI_twice)*M_PI_twice - M_PI;		// wrap p to -pi..pi		
//			cp = ( p<M_PI_half )  ? p + M_PI_half : p - M_PI_3half;

			//pre-fetch next data
			csx = cs.c1[j+1];
			csy = cs.c2[j+1];
			csz = cs.c3[j+1];	

//			Ar += sfs[j]*sign_sin*sine(p);
//			Ai += sfs[j]*sine(cp);
			Ar += sfs[j]*cos(p);
			Ai += sfs[j]*sin(p);			
		}
		a(localframe,0)=complex<double>(Ar,Ai); // sum at this location
}

void AllScatterDevice::scatter_frames_norm1(CartesianCoor3D& q) {
	for(size_t i = 0; i < a.size1(); ++i)
	{
		timer.start("sd:fs:f");												
		scatter_frame_norm1(i,q);
		timer.stop("sd:fs:f");														
	}
}


std::vector<std::complex<double> > AllScatterDevice::correlate_frames_fftw() {
	// each node has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

	// negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();

	timer.start("sd:corr::areduce");										
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());
	timer.stop("sd:corr::areduce");										

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		local_A[ci] = a(ci,0);
	}

	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; 
	if (p_thisworldcomm->rank()==0) {
    	all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());		    
	}

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:corr:gather");										
	boost::mpi::gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0],0);
	timer.stop("sd:corr:gather");										

    if (p_thisworldcomm->rank()==0) {

    	EvenDecompose edecomp(p_sample->coordinate_sets.size(),p_thisworldcomm->size());

    	// this has interleaving A , have to be indexed away
    	vector<complex<double> > A; A.resize(p_sample->coordinate_sets.size());	

    	for(size_t i = 0; i < p_thisworldcomm->size(); ++i)
    	{
    		vector<size_t> findexes = edecomp.indexes_for(i);
    		for(size_t j = 0; j < findexes.size(); ++j)
    		{
    			A[ findexes[j] ] = complex<double>(all_Ar[ 2*maxCSsize*i + 2*j ],all_Ar[ 2*maxCSsize*i + 2*j +1]);
    		}
    	}

        // A contains all complex values in the right order.
        
    	vector<complex<double> > correlated_a; correlated_a.resize(p_sample->coordinate_sets.size(),complex<double>(0,0));        
    	timer.start("sd:cor:correlate");
        Complex *f1=FFTWComplex(p_sample->coordinate_sets.size()*2);
        Complex *f2=FFTWComplex(p_sample->coordinate_sets.size()*2);

        for(size_t i = 0; i < p_sample->coordinate_sets.size(); ++i) {
           f1[i]=A[i];
        } 
        for(size_t i = p_sample->coordinate_sets.size(); i < p_sample->coordinate_sets.size()*2; ++i) {
           f1[i]=0;
        } 
        

        fft1d Forward(p_sample->coordinate_sets.size()*2,-1);
        
        fft1d Backward(p_sample->coordinate_sets.size()*2,1);

        Forward.fft(f1);
        for(size_t i = 0; i < p_sample->coordinate_sets.size()*2; ++i) {
           f2[i]=f1[i]*conj(f1[i]); 
        } 
        Backward.fft(f2);

        double cssize = p_sample->coordinate_sets.size();

        // normalize by:
        // 2 * cssize = FFTW introduced scaling 
        // cssize - i = compensate sampling weight (like in direct correlation)
        for(size_t i = 0; i < p_sample->coordinate_sets.size(); ++i) {           
           correlated_a[i]=f2[i]*( 1.0 / ( 2*cssize * (cssize -i ) ) ); 
        }
         
        FFTWdelete(f1);
        FFTWdelete(f2);              
  
    	timer.stop("sd:corr:correlate");
        
        return correlated_a;
    } else {
        vector<std::complex<double> > out; // return an empty array for none-root
        return out;
    }
}


std::vector<std::complex<double> > AllScatterDevice::correlate_frames() {
	// each nodes has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

    size_t NF = p_sample->coordinate_sets.size();
    size_t NN = p_thisworldcomm->size();
    
	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize = CSsize;
	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:corr::areduce");										
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());
	timer.stop("sd:corr::areduce");										

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < CSsize; ++ci)
	{
		local_A[ci] = a(ci,0);
	}


	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*NN);	

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:corr:agather");										
	boost::mpi::all_gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0]);
	timer.stop("sd:corr:agather");										

	EvenDecompose edecomp(NF,NN);

	// this has interleaving A , have to be indexed away
	vector<complex<double> > A; A.resize(NF);	

	for(size_t i = 0; i < p_thisworldcomm->size(); ++i)
	{
		vector<size_t> findexes = edecomp.indexes_for(i);
		for(size_t j = 0; j < findexes.size(); ++j)
		{
			A[ findexes[j] ] = complex<double>(all_Ar[ 2*maxCSsize*i + 2*j ],all_Ar[ 2*maxCSsize*i + 2*j +1]);
		}
	}

	RModuloDecompose rmdecomp(NF,NN);
	vector<size_t> mysteps = rmdecomp.indexes_for(p_thisworldcomm->rank());
	
	vector<complex<double> > correlated_a; correlated_a.resize(NF,complex<double>(0,0));

	timer.start("sd:cor:correlate");
	
    complex<double> mean; mean =0.0;
	if (Params::Inst()->scattering.correlation.zeromean) {
	    for(size_t i = 0; i < NF; ++i)
	    {
            mean += A[i];
	    }
        mean = mean * (1.0/NF);
	}
												
	for(size_t i = 0; i < mysteps.size(); ++i)
	{
		size_t tau = mysteps[i];
		size_t last_starting_frame = NF-tau;
		for(size_t k = 0; k < last_starting_frame; ++k) // this iterates the starting frame
		{
			complex<double> a1 = A[k] - mean;
			complex<double> a2 = A[k+tau] - mean;
			correlated_a[tau] += conj(a1)*a2;
		}
		correlated_a[tau] /= (last_starting_frame); // maybe a conj. multiply here		
	}
	timer.stop("sd:corr:correlate");											

	vector<double> ain_r = flatten(correlated_a);
	
	vector<double> aout_r; aout_r.resize(ain_r.size());

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:corr:reduce");											
	boost::mpi::reduce(*p_thisworldcomm,&ain_r[0],ain_r.size(),&aout_r[0],std::plus<double>(),0);
	timer.stop("sd:corr:reduce");											

	return compress(aout_r);
}

vector<complex<double> > AllScatterDevice::conjmultiply_frames() {
	// each node has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)


	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;
	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:conmul:areduce");										
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());
	timer.stop("sd:conmul:areduce");										

	// in conjmultiply we just have to conj multiply each a(x,0)
	complex<double> cA = 0; 
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		a(ci,0) = conj(a(ci,0))*a(ci,0);
	}

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		local_A[ci] = a(ci,0);
	}

	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());	

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:conjmul:gather");										
	boost::mpi::gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0],0);
	timer.stop("sd:conjmul:gather");										

	if (p_thisworldcomm->rank()==0) {
		EvenDecompose edecomp(p_sample->coordinate_sets.size(),p_thisworldcomm->size());

		// this has interleaving A , have to be indexed away
		vector<complex<double> > A; A.resize(p_sample->coordinate_sets.size());	

		for(size_t i = 0; i < p_thisworldcomm->size(); ++i)
		{
			vector<size_t> findexes = edecomp.indexes_for(i);
			for(size_t j = 0; j < findexes.size(); ++j)
			{
				A[ findexes[j] ] = complex<double>(all_Ar[ 2*maxCSsize*i + 2*j ],all_Ar[ 2*maxCSsize*i + 2*j +1]);
			}
		}

		return A;
	} else {
		// return is empty
		vector<complex<double> > A;
		return A;
	}
	
}

void AllScatterDevice::superpose_spectrum(vector<complex<double> >& spectrum, vector<complex<double> >& fullspectrum) {
	for(size_t j = 0; j < spectrum.size(); ++j)
	{
		fullspectrum[j] += spectrum[j];
	}
}


void AllScatterDevice::execute(CartesianCoor3D q) {
			
	vector<CartesianCoor3D> qvectors;	
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
		if (Params::Inst()->scattering.average.orientation.vectors.type=="sphere") {
			double ql = q.length();
			
			for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
			{
				qvectors.push_back(ql*Params::Inst()->scattering.average.orientation.vectors[i]);
			}					
		} else if (Params::Inst()->scattering.average.orientation.vectors.type=="cylinder") {
			CartesianCoor3D o = Params::Inst()->scattering.average.orientation.vectors.axis;
			// make sure o is normalized;
			o = o / o.length();
			
			// get the part of the scattering vector perpenticular to the o- orientation
			CartesianCoor3D qparallel = (o*q)*o; 
			CartesianCoor3D qperpenticular = q - qparallel; 			
			double qperpenticular_l = qperpenticular.length();
			double qparallel_l = qparallel.length();

			CartesianCoor3D e1 = o.cross_product(qperpenticular) ;
			CartesianCoor3D e2 = qperpenticular;
			
			if (qperpenticular_l==0.0)  { 
				qvectors.push_back( qparallel );
			}
			else {
				for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
				{
					CartesianCoor3D& vec = Params::Inst()->scattering.average.orientation.vectors[i];
					CartesianCoor3D qnew = qparallel + vec.x * e1 + vec.y * e2;					
					qvectors.push_back(qnew);
				}
			}				
		}
	} else {
		qvectors.push_back(q);
	}
				
	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; spectrum.resize(p_sample->coordinate_sets.size());

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	for(size_t qi = 0; qi < qvectors.size(); ++qi) // blocking of N
	{
		timer.start("sd:fs");
		scatter_frames_norm1(qvectors[qi]); // put summed scattering amplitudes into first atom entry		
		timer.stop("sd:fs");
		
		vector<complex<double> > thisspectrum;
		if (Params::Inst()->scattering.correlation.type=="time") {
			if (Params::Inst()->scattering.correlation.method=="direct") {
				timer.start("sd:correlate");					
				thisspectrum = correlate_frames(); // if correlation, otherwise do a elementwise conj multiply here			
				timer.stop("sd:correlate");					
    		} else if (Params::Inst()->scattering.correlation.method=="fftw") {
    				timer.start("sd:correlate");
    				thisspectrum = correlate_frames_fftw(); // if correlation, otherwise do a elementwise conj multiply here			
    				timer.stop("sd:correlate");					
    		} else {
				Err::Inst()->write("Correlation method not understood. Supported methods: direct fftw");
				throw;
			}
		} else {
			timer.start("sd:conjmul");				
			thisspectrum = conjmultiply_frames();
			timer.stop("sd:conjmul");									
		}
		
		if (p_thisworldcomm->rank()==0) {
			timer.start("sd:superpose");										
			superpose_spectrum(thisspectrum,spectrum);
			timer.stop("sd:superpose");														
		}											
	}
	
	
	for(size_t si = 0; si < spectrum.size(); ++si)
	{
		spectrum[si] /= qvectors.size() ;
	}
	
	m_spectrum = spectrum;

}

vector<complex<double> >& AllScatterDevice::get_spectrum() {
	return m_spectrum;
}
