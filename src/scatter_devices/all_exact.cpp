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
#include <boost/math/special_functions.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "control.hpp"
#include "sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"
#include "vector_unfold.hpp"

using namespace std;

AllExactScatterDevice::AllExactScatterDevice(boost::mpi::communicator& thisworld, Sample& sample) {

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t NN = thisworld.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	size_t rank = thisworld.rank();
	EvenDecompose edecomp(NF,NN);
	
	myframes = edecomp.indexes_for(rank);
	
	if (thisworld.rank()==0) {	

        size_t memusage_scatmat = 2*sizeof(double)*myframes.size()*1;
                
        size_t memusage_per_cs = 3*sizeof(double)*NA;
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
	
	p_sample->coordinate_sets.set_representation(SPHERICAL);		
	p_sample->coordinate_sets.set_selection(p_sample->atoms.selections[target]);
	if (Params::Inst()->scattering.center) {
    	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	}
	
	scatterfactors.set_sample(*p_sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

// acts like scatter_frame, but does the summation in place
void AllExactScatterDevice::scatter_frame_norm1(size_t iframe, CartesianCoor3D& q) {

	size_t noa = p_sample->coordinate_sets.get_selection().indexes.size();
	
	timer.start("sd:fs:f:ld");	
	CoordinateSet& cs = p_sample->coordinate_sets.load(myframes[iframe]); 
	timer.stop("sd:fs:f:ld");	
	vector<double>& sfs = scatterfactors.get_all();

	using namespace boost::numeric::ublas::detail;
	
	double ql = q.length();
	
    double A=0;
    for(size_t i = 0; i < noa; ++i)
    {
		double r1   = cs.c1[i];
        double esf1 = sfs[i];
        A += 2*esf1*esf1;
	    for(size_t j = i+1; j < noa; ++j) {
    		double r2   = cs.c1[j];
            double rd = r1-r2;
            // *2 b/c of symmetry
            A += 2*esf1*sfs[j]*boost::math::sinc_pi(ql*rd);
	    }
    }
  
	a(iframe,0)=complex<double>(A,0); // sum at this location
}


void AllExactScatterDevice::scatter_frames_norm1(CartesianCoor3D& q) {
	
	for(size_t i = 0; i < myframes.size(); ++i)
	{
		timer.start("sd:fs:f");			    
		scatter_frame_norm1(i,q);
		timer.stop("sd:fs:f");				
	}
}

vector<complex<double> > AllExactScatterDevice::gather_frames() {
	// each node has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;
    timer.start("sd:gf:areduce");	
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());
    timer.stop("sd:gf:areduce");

	// for multipole we already have the conj-multiplied version 

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		local_A[ci] = a(ci,0);
	}

	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());	

    timer.start("sd:gf:gather");	
	boost::mpi::gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0],0);
    timer.stop("sd:gf:gather");	

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
		vector<complex<double> > A;
		return A;
	}
}

void AllExactScatterDevice::superpose_spectrum(vector<complex<double> >& spectrum, vector<complex<double> >& fullspectrum) {
	for(size_t j = 0; j < spectrum.size(); ++j)
	{
		fullspectrum[j] += spectrum[j];
	}
}

void AllExactScatterDevice::multiply_alignmentfactors(CartesianCoor3D q) {
    // not yet implemented
    Err::Inst()->write("centering w/ alignment not yet implemented");
    throw;
}

void AllExactScatterDevice::execute(CartesianCoor3D q) {
			
	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; spectrum.resize(p_sample->coordinate_sets.size());

	timer.start("sd:sf:update");
	scatterfactors.update(q); 
	timer.stop("sd:sf:update");
	
	timer.start("sd:fs");	    
	scatter_frames_norm1(q); // put summed scattering amplitudes into first atom entry
	timer.stop("sd:fs");
	
	if (Params::Inst()->scattering.center) {
        multiply_alignmentfactors(q);
	}	
		
	if (Params::Inst()->scattering.correlation.type=="time") {
		Err::Inst()->write("Correlation not supported with the exact method for spherical averaging");
		throw;
	} else {
		timer.start("sd:gatherframes");			    
		spectrum = gather_frames();
		timer.stop("sd:gatherframes");				
	}

	m_spectrum = spectrum;

}

vector<complex<double> >& AllExactScatterDevice::get_spectrum() {
	return m_spectrum;
}


