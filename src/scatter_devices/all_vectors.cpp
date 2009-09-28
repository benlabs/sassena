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
#include "log.hpp"
#include "parameters.hpp"
#include "database.hpp"
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
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
}

// acts like scatter_frame, but does the summation in place
void AllScatterDevice::scatter_frame_norm1(size_t localframe, CartesianCoor3D& q) {

	size_t noa = coordinate_sets.get_selection().size();
	
	timer.start("sd:fs:f:ld");												
	CoordinateSet& cs = coordinate_sets.load(myframes[localframe]);
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

		double csx = cs.x[0];
		double csy = cs.y[0];
		double csz = cs.z[0];


		for(size_t j = 0; j < noa; ++j) {
			p =  csx*q.x + csy*q.y + csz*q.z;
//			double sign_sin = (p<0) ? -1.0 : 1.0;
//			ap = abs(p);
//			p = ap - long(ap/M_PI_twice)*M_PI_twice - M_PI;		// wrap p to -pi..pi		
//			cp = ( p<M_PI_half )  ? p + M_PI_half : p - M_PI_3half;

			//pre-fetch next data
			csx = cs.x[j+1];
			csy = cs.y[j+1];
			csz = cs.z[j+1];	

//			Ar += sfs[j]*sign_sin*sine(p);
//			Ai += sfs[j]*sine(cp);
			Ar += sfs[j]*sin(p);
			Ai += sfs[j]*cos(p);			
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

std::vector<std::complex<double> > AllScatterDevice::correlate_frames() {
	// each nodes has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

	//negotiate maximum size for coordinatesets
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
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());	

	if (Params::Inst()->debug.barriers) p_thisworldcomm->barrier();
	timer.start("sd:corr:agather");										
	boost::mpi::all_gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0]);
	timer.stop("sd:corr:agather");										

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

	RModuloDecompose rmdecomp(p_sample->frames.size(),p_thisworldcomm->size());
	vector<size_t> mysteps = rmdecomp.indexes_for(p_thisworldcomm->rank());
	
	vector<complex<double> > correlated_a; correlated_a.resize(p_sample->frames.size(),complex<double>(0,0));

	timer.start("sd:cor:correlate");											
	for(size_t i = 0; i < mysteps.size(); ++i)
	{
		size_t tau = mysteps[i];
		size_t last_starting_frame = p_sample->frames.size()-tau;
		for(size_t k = 0; k < last_starting_frame; ++k) // this iterates the starting frame
		{
			complex<double> a1 = A[k];
			complex<double> a2 = A[k+tau];
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

void AllScatterDevice::execute(CartesianCoor3D& q) {
			
	string avm = Params::Inst()->scattering.average.orientation.method;
				
	timer.start("sd:vector:unfold");
	VectorUnfold* p_vectorunfold = NULL;			
	if (avm=="bruteforce") {
		string avv = Params::Inst()->scattering.average.orientation.vectors;
		string avt = Params::Inst()->scattering.average.orientation.type;

		if (avt=="sphere") {
			SphereVectorUnfold* p_vu = new SphereVectorUnfold(q);
			p_vu->set_resolution(Params::Inst()->scattering.average.orientation.resolution);
			p_vu->set_vectors(Params::Inst()->scattering.average.orientation.vectors);
			p_vu->set_seed(0);
			p_vectorunfold = p_vu; 
		}
		if (avt=="cylinder") {
			CylinderVectorUnfold* p_vu = new CylinderVectorUnfold(q);
			p_vu->set_resolution(Params::Inst()->scattering.average.orientation.resolution);
			p_vu->set_vectors(Params::Inst()->scattering.average.orientation.vectors);
			p_vu->set_axis(Params::Inst()->scattering.average.orientation.axis);
			p_vu->set_seed(0);			
			p_vectorunfold = p_vu; 
		}		 
		else if (avt=="file") {
			p_vectorunfold = new FileVectorUnfold(q);
		}
	}
	else if (avm=="multipole") {
		p_vectorunfold = new NoVectorUnfold(q);
	}
	else { // default: no unfold
		p_vectorunfold = new NoVectorUnfold(q);
	}

	p_vectorunfold->execute();
	timer.stop("sd:vector:unfold");

	vector<CartesianCoor3D>& qvectors = p_vectorunfold->vectors();

	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; spectrum.resize(p_sample->frames.size());

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	for(size_t qi = 0; qi < qvectors.size(); ++qi)
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
			} else {
				Err::Inst()->write("Correlation method not understood. Supported methods: direct");
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
	
	delete p_vectorunfold;
}

vector<complex<double> >& AllScatterDevice::get_spectrum() {
	return m_spectrum;
}
