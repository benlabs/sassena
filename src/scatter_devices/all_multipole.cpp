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
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"
#include "sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"
#include "vector_unfold.hpp"

using namespace std;

AllMScatterDevice::AllMScatterDevice(boost::mpi::communicator& thisworld, Sample& sample) {

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t nn = thisworld.size(); // Number of Nodes
	size_t na = sample.atoms.selections[target].size(); // Number of Atoms
	size_t nf = sample.frames.size();

	size_t rank = thisworld.rank();

	EvenDecompose edecomp(nf,nn);
	
	myframes = edecomp.indexes_for(rank);
	if (thisworld.rank()==0) {	
		Info::Inst()->write(string("Memory Requirements per node: "));
		Info::Inst()->write(string("Scattering Matrix: ")+to_s(2*8*myframes.size()*1)+string(" bytes"));

		size_t numframes = myframes.size(); 
		Info::Inst()->write(string("Coordinate Sets: ")+to_s(3*8*myframes.size()*na)+string(" bytes"));


		if (Params::Inst()->limits.coordinatesets_cache_max!=0 && Params::Inst()->limits.coordinatesets_cache_max<numframes) {
			Warn::Inst()->write(string("Insufficient Buffer size for coordinate sets. This is a huge bottleneck for performance!"));
			Warn::Inst()->write(string("Need at least: ")+to_s(numframes));
		}		
	}
	a.resize(myframes.size(),1); // n x 1 = frames x 1
	
	
	
	coordinate_sets.set_sample(sample);
	coordinate_sets.set_selection(sample.atoms.selections[target]);
	sample.atoms.assert_selection(Params::Inst()->scattering.average.orientation.origin);
	coordinate_sets.set_origin(sample.atoms.selections[Params::Inst()->scattering.average.orientation.origin]);
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

// acts like scatter_frame, but does the summation in place
void AllMScatterDevice::scatter_frame_norm1(size_t iframe, CartesianCoor3D& q) {

	size_t noa = coordinate_sets.get_selection().size();
	
	// this is a specially hacked version of CoordinateSet , it contains r, phi, theta at x,y,z repectively
	CoordinateSet& cs = coordinate_sets.load(iframe); 
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



	using namespace boost::numeric::ublas::detail;
	
	long resolution = Params::Inst()->scattering.average.orientation.resolution;	
	const int lmax=int(resolution);
	
	std::vector<std::vector<complex<double> > > almv;

  	almv.resize(lmax+1);
  	for (int l=0;l<=lmax;l++) {
   		almv[l].resize(2*l+1,complex<double>(0,0));
  	}

	double ql = q.length();
	double M_PI_four = 4*M_PI;

	for(size_t j = 0; j < noa; ++j) {
		// again... this is a special hack: r=x,phi=y,theta=z
		double r   = cs.x[j];
	 	double phi = cs.y[j];
		double theta = cs.z[j];
		
		double esf = sfs[j];
		
		for (int l=0;l<=lmax;l++) {
			complex<double> fmpiilesf = M_PI_four*pow(complex<double>(0,1.0),l) * esf;
			double aabess = boost::math::sph_bessel(l,ql*r);
		
			for (int m=-l;m<=l;m++) {
		
			complex<double> aa = conj(boost::math::spherical_harmonic(l,m,theta,phi)); 

			almv[l][m+l] += fmpiilesf * aabess* aa;
			}
		}
	}

	// we need to multiply w/ (lmax+1)^2 here b/c single ampltiudes are average-summed
//	double lmax1sqr = (lmax+1) / sqrt(M_PI_four); // not necessary here.....
	complex<double> A(0,0);
	for (int l=0;l<=lmax;l++) {
		for (int m=-l;m<=l;m++) { 
//			A+= lmax1sqr*(almv[l][m+l]) * lmax1sqr*conj(almv[l][m+l]) ; 
			A+= (almv[l][m+l]) * conj(almv[l][m+l]) ; 
		}
	}
	
	A /= M_PI_four;
	a(iframe,0)=A; // sum at this location
}

void AllMScatterDevice::scatter_frames_norm1(CartesianCoor3D& q) {
	
	for(size_t i = 0; i < myframes.size(); ++i)
	{
		scatter_frame_norm1(i,q);
	}
}

vector<complex<double> > AllMScatterDevice::gather_frames() {
	// each node has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());


	// for multipole we already have the conj-multiplied version 

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		local_A[ci] = a(ci,0);
	}

	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());	

	boost::mpi::gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0],0);

	if (p_thisworldcomm->rank()==0) {
		EvenDecompose edecomp(p_sample->frames.size(),p_thisworldcomm->size());

		// this has interleaving A , have to be indexed away
		vector<complex<double> > A; A.resize(p_sample->frames.size());	

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
	}
}

void AllMScatterDevice::superpose_spectrum(vector<complex<double> >& spectrum, vector<complex<double> >& fullspectrum) {
	for(size_t j = 0; j < spectrum.size(); ++j)
	{
		fullspectrum[j] += spectrum[j];
	}
}

void AllMScatterDevice::execute(CartesianCoor3D& q) {
			
	string avm = Params::Inst()->scattering.average.orientation.method;
				
	VectorUnfold* p_vectorunfold = NULL;							
	p_vectorunfold = new NoVectorUnfold(q);	
	p_vectorunfold->execute();
	vector<CartesianCoor3D>& qvectors = p_vectorunfold->vectors();

	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; spectrum.resize(p_sample->frames.size());

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	for(size_t qi = 0; qi < qvectors.size(); ++qi)
	{		
		scatter_frames_norm1(qvectors[qi]); // put summed scattering amplitudes into first atom entry
	
		vector<complex<double> > thisspectrum;
		if (Params::Inst()->scattering.correlation.type=="time") {
			Err::Inst()->write("Correlation not supported with the multipole method for spherical averaging");
			throw;
		} else {
			thisspectrum = gather_frames();
		}

		if (p_thisworldcomm->rank()==0) superpose_spectrum(thisspectrum,spectrum);
	}
	
	for(size_t si = 0; si < spectrum.size(); ++si)
	{
		spectrum[si] /= qvectors.size();
	}
	
	m_spectrum = spectrum;
	
	delete p_vectorunfold;
}

vector<complex<double> >& AllMScatterDevice::get_spectrum() {
	return m_spectrum;
}
