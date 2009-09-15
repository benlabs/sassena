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

vector<double> flatten(vector<complex<double> >& cvalues) {
	vector<double> rvalues(cvalues.size()*2);
	for(size_t i = 0; i < cvalues.size(); ++i)
	{
		rvalues[2*i]=cvalues[i].real();
		rvalues[2*i+1]=cvalues[i].imag();
	}
	return rvalues;
}

vector<complex<double> > compress(vector<double>& rvalues) {
	vector<complex<double> > cvalues(rvalues.size()/2);
	for(size_t i = 0; i < cvalues.size(); ++i)
	{
		cvalues[i]=complex<double>(rvalues[2*i],rvalues[2*i+1]);
	}
	return cvalues;
}

SelfScatterDevice::SelfScatterDevice(boost::mpi::communicator& thisworld, Sample& sample){

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;
	

	size_t nn = thisworld.size(); // Number of Nodes
	size_t na = sample.atoms.selections[target].size(); // Number of Atoms
	size_t nf = sample.frames.size();

	size_t rank = thisworld.rank();

	EvenDecompose edecomp(nf,nn);

	if (nn>na) {
		Err::Inst()->write("Not enough atoms for decomposition");
		throw;
	}

	size_t minatoms = na/nn;	
	size_t leftoveratoms = na % nn;	

//	CoordinateSet* p_cs =NULL;

	size_t* p_indexes; size_t* p_myindexes;
	// construct particle_trajectories first
	vector<size_t> selectionindexes(sample.atoms.selections[target].size());
	if (rank==0) {
//		sample.frames.load(0,sample.atoms);
//		p_cs = new CoordinateSet(sample.frames.current(),sample.atoms.selections[target]);
		// don't use the real atom indexes, but the selection indexes instead, this will work well with the scatterfactors class			
//		p_indexes = &(sample.atoms.selections[target][0]); 
		for(size_t i = 0; i < selectionindexes.size(); ++i) selectionindexes[i]=i;
		p_indexes = &(selectionindexes[0]);
	}

	vector<size_t> myindexes(minatoms); p_myindexes = &(myindexes[0]);
	
	// this is an atom decomposition!
	boost::mpi::scatter(thisworld,p_indexes,p_myindexes,minatoms,0);	
	for(size_t i = 0; i < minatoms; ++i)
	{
		particle_trajectories.push_back(ParticleTrajectory(myindexes[i]) );
	}

	if (leftoveratoms>0) {
		
		if (rank==0) {
			p_indexes = &(sample.atoms.selections[target][minatoms]);
		}

		size_t lastindex; p_myindexes = &(lastindex);
		boost::mpi::scatter(thisworld,p_indexes,p_myindexes,1,0);	
		if (rank<leftoveratoms) {
			particle_trajectories.push_back(ParticleTrajectory(lastindex) );		
		}	
	}

	for(size_t r = 0; r < nn; ++r)
	{

		vector<size_t> assigned_frames = edecomp.indexes_for(r);
		
		for(size_t i = 0; i < assigned_frames.size(); ++i)
		{
			double  *p_xc,*p_yc,*p_zc;	
			
			CoordinateSet* p_cs =NULL;
			if (rank==r) {
				sample.frames.load(assigned_frames[i],sample.atoms);
				p_cs = new CoordinateSet(sample.frames.current(),sample.atoms.selections[target]);
				
				p_xc = &(p_cs->x[0]);
				p_yc = &(p_cs->y[0]);
				p_zc = &(p_cs->z[0]);
							
			}
			
			vector<double> myxc(minatoms);
			vector<double> myyc(minatoms);
			vector<double> myzc(minatoms);

			boost::mpi::scatter(thisworld,p_xc,&(myxc[0]),minatoms,r);
			boost::mpi::scatter(thisworld,p_yc,&(myyc[0]),minatoms,r);
			boost::mpi::scatter(thisworld,p_zc,&(myzc[0]),minatoms,r);	

			for(size_t pi = 0; pi < minatoms; ++pi)
			{
				particle_trajectories[pi].push_back(CartesianCoor3D(myxc[pi],myyc[pi],myzc[pi]));
			}		

			if (rank==r) {
				p_xc = &(p_cs->x[minatoms]);
				p_yc = &(p_cs->y[minatoms]);
				p_zc = &(p_cs->z[minatoms]);
			}

			double mylastx; 
			double mylasty; 
			double mylastz; 
			
			boost::mpi::scatter(thisworld,p_xc,&mylastx,1,r);
			boost::mpi::scatter(thisworld,p_yc,&mylasty,1,r);
			boost::mpi::scatter(thisworld,p_zc,&mylastz,1,r);	

			if (rank<leftoveratoms) {
				particle_trajectories[minatoms].push_back(CartesianCoor3D(mylastx,mylasty,mylastz));
			}
				
			if (rank==r) delete p_cs;

		}
	}
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
	a.resize(particle_trajectories.size(),sample.frames.size());

}

void SelfScatterDevice::scatter_particle(size_t iparticle, CartesianCoor3D& q) {
	size_t row = iparticle;
	ParticleTrajectory& thisparticle = particle_trajectories[iparticle];
	// this is broken
	double s = scatterfactors.get(thisparticle.atom_selection_index());
		
	for(size_t j = 0; j < a.size2(); ++j)
	{
		double x1 = thisparticle.x[j];
		double y1 = thisparticle.y[j];
		double z1 = thisparticle.z[j];	
		double p1 = x1*q.x+y1*q.y+z1*q.z;

		a(iparticle,j) = exp(-1.0*complex<double>(0,p1))*s;
	}
	
//	double Ar = 0.0;
//	double Ai = 0.0;
//	double p, cp, ap;
//	double M_PI_half = M_PI/2;
//	double M_PI_3half = 3*M_PI/2;	
//	double M_PI_twice = 2*M_PI;
//
//	double px = thisparticle.x[0];
//	double py = thisparticle.y[0];
//	double pz = thisparticle.z[0];
//
//	for(size_t j = 0; j < a.size2(); ++j) {
//		p =  px*q.x + py*q.y + pz*q.z;
//		double sign_sin = (p<0) ? -1.0 : 1.0;
//		ap = abs(p);
//		p = ap - long(ap/M_PI_twice)*M_PI_twice - M_PI;		// wrap p to -pi..pi		
//		cp = ( p<M_PI_half )  ? p + M_PI_half : p - M_PI_3half;
//
//		//pre-fetch next data
//		px = thisparticle.x[j+1];
//		py = thisparticle.y[j+1];
//		pz = thisparticle.z[j+1];	
//
//		Ar = sfs[j]*sign_sin*sine(p);
//		Ai = sfs[j]*sine(cp);
//	}
//	a(iframe,0)=complex<double>(Ar,Ai); // sum at this location	
	
}

void SelfScatterDevice::scatter_particles(CartesianCoor3D& q) {
	for(size_t i = 0; i < a.size1(); ++i)
	{
		scatter_particle(i,q);
	}
}

void SelfScatterDevice::scatter_particle_instant(size_t iparticle, CartesianCoor3D& q) {
	ParticleTrajectory& thisparticle = particle_trajectories[iparticle];
	// this is broken
	double s = scatterfactors.get(thisparticle.atom_selection_index());
		
	for(size_t j = 0; j < a.size2(); ++j)
	{
		a(iparticle,j) = s;
	}	
}

void SelfScatterDevice::scatter_particles_instant(CartesianCoor3D& q) {
	for(size_t i = 0; i < a.size1(); ++i)
	{
		scatter_particle_instant(i,q);
	}
}

void SelfScatterDevice::correlate_particles() {

	for(size_t i = 0; i < a.size1(); ++i)
	{
		for(size_t j = 0; j < a.size2(); ++j) // this iterates tau
		{
			complex<double> a1conj = conj( a (i,j) );
			complex<double> asum(0,0);
			size_t last_starting_frame = a.size2()-j;
			for(size_t k = 0; k < last_starting_frame; ++k) // this iterates the starting frame
			{
				asum += a1conj*( a (i,k+j) ) ;
			}
			a (i,j)  = asum;
			a (i,j) /= last_starting_frame; // store it in (i,j)
		}
		
	}
	// now a contains correlated values
}

void SelfScatterDevice::conjmultiply_particles() {

	for(size_t i = 0; i < a.size1(); ++i)
	{		
		for(size_t j = 0; j < a.size2(); ++j) // this iterates t
		{			
			a(i,j) *= conj( a (i,j) );
		}
	}
	// now a contains conjmultiplied values
}

void SelfScatterDevice::crosssum_particles() {
	for(size_t i = 1; i < a.size1(); ++i)
	{
		for(size_t j = 0; j < a.size2(); ++j)
		{
			a(0,j) += a(i,j);
		}
	}
	// now a(0,x) contains the summed spectrum
}

void SelfScatterDevice::execute(CartesianCoor3D& q) {
	
   string avm = Params::Inst()->scattering.average.method;
   
   VectorUnfold* p_vectorunfold = NULL;			
   if (avm=="bruteforce") {
   	string avv = Params::Inst()->scattering.average.vectors;
   	string avt = Params::Inst()->scattering.average.type;
   
   	if (avt=="sphere") {
   		SphereVectorUnfold* p_vu = new SphereVectorUnfold(q);
   		p_vu->set_resolution(Params::Inst()->scattering.average.resolution);
   		p_vu->set_vectors(Params::Inst()->scattering.average.vectors);
   		p_vu->set_seed(0);
   		p_vectorunfold = p_vu; 
   	}
   	if (avt=="cylinder") {
   		CylinderVectorUnfold* p_vu = new CylinderVectorUnfold(q);
   		p_vu->set_resolution(Params::Inst()->scattering.average.resolution);
   		p_vu->set_vectors(Params::Inst()->scattering.average.vectors);
   		p_vu->set_axis(Params::Inst()->scattering.average.axis);
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
   
   vector<CartesianCoor3D>& qvectors = p_vectorunfold->qvectors();
	
	
	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; 
	
	if (p_thisworldcomm->rank()==0)  spectrum.resize(a.size2()); // only rank0 node aggregrates results

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	for(size_t qi = 0; qi < qvectors.size(); ++qi)
	{
		vector<complex<double> > thisspectrum;
		if (Params::Inst()->scattering.correlation.type=="time") {
			scatter_particles(qvectors[qi]);

			if (Params::Inst()->scattering.correlation.method=="direct") {
				correlate_particles(); // if correlation, otherwise do a elementwise conj multiply here	
			} else {
				Err::Inst()->write("Correlation method not understood. Supported methods: direct");
				throw;
			}
		} else {
			scatter_particles_instant(qvectors[qi]); // this version takes advantage of a*a(conj) = (real,0)
			conjmultiply_particles();
		}
		crosssum_particles(); // the spectrum is in a(0,x)		
		assemble_spectrum(); // now rank==0 has the full spectrum in a(0,x)


		if (p_thisworldcomm->rank()==0) superpose_spectrum(spectrum);
	}
	
	for(size_t si = 0; si < spectrum.size(); ++si)
	{
		spectrum[si] /= qvectors.size();
	}

	m_spectrum = spectrum;
}

vector<complex<double> >& SelfScatterDevice::get_spectrum() {
	return m_spectrum;
}

void SelfScatterDevice::assemble_spectrum() {

	// local spectrum: a(0,x)
	
	// aggregate everything on rank 0
	size_t rank = p_thisworldcomm->rank();
	size_t nn = p_thisworldcomm->size();
	
	a.size2();
	complex<double>* pAin = &a(0,0);
	vector<double> Ain;
	vector<double> Aout; Aout.resize(a.size2()*2,0);
	for(size_t i = 0; i < a.size2(); ++i)
	{
		Ain.push_back(a(0,i).real());
		Ain.push_back(a(0,i).imag());		
	}
	
	boost::mpi::reduce(*p_thisworldcomm,&Ain[0],2*a.size2(),&Aout[0],std::plus<double>(),0);
	
	if (p_thisworldcomm->rank()==0) {
		for(size_t i = 0; i < a.size2(); ++i)
		{
			a(0,i) = complex<double>(Aout[2*i],Aout[2*i+1]);		
		}
	}
}

void SelfScatterDevice::superpose_spectrum(vector<complex<double> >& fullspectrum) {
	for(size_t j = 0; j < a.size2(); ++j)
	{
		fullspectrum[j] += a(0,j);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

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
	
	a.resize(myframes.size(),sample.atoms.selections[target].size()); // n x m = frames x atoms
	
	coordinate_sets.set_sample(sample);
	coordinate_sets.set_selection(sample.atoms.selections[target]);
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

// acts like scatter_frame, but does the summation in place
void AllScatterDevice::scatter_frame_norm1(size_t iframe, CartesianCoor3D& q) {
	
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



		double Ar = 0.0;
		double Ai = 0.0;
		double p, cp, ap;
		double M_PI_half = M_PI/2;
		double M_PI_3half = 3*M_PI/2;	
		double M_PI_twice = 2*M_PI;

		double csx = cs.x[0];
		double csy = cs.y[0];
		double csz = cs.z[0];

		for(size_t j = 0; j < a.size2(); ++j) {
			p =  csx*q.x + csy*q.y + csz*q.z;
			double sign_sin = (p<0) ? -1.0 : 1.0;
			ap = abs(p);
			p = ap - long(ap/M_PI_twice)*M_PI_twice - M_PI;		// wrap p to -pi..pi		
			cp = ( p<M_PI_half )  ? p + M_PI_half : p - M_PI_3half;

			//pre-fetch next data
			csx = cs.x[j+1];
			csy = cs.y[j+1];
			csz = cs.z[j+1];	

			Ar += sfs[j]*sign_sin*sine(p);
			Ai += sfs[j]*sine(cp);
		}
		a(iframe,0)=complex<double>(Ar,Ai); // sum at this location
}

void AllScatterDevice::scatter_frames_norm1(CartesianCoor3D& q) {
	
	for(size_t i = 0; i < myframes.size(); ++i)
	{
		scatter_frame_norm1(i,q);
	}
}

void AllScatterDevice::scatter_frame(size_t iframe, CartesianCoor3D& q) {

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
	
	
	
	double Ar = 0.0;
	double Ai = 0.0;
	double p, cp, ap;
	double M_PI_half = M_PI/2;
	double M_PI_3half = 3*M_PI/2;	
	double M_PI_twice = 2*M_PI;
	
	double csx = cs.x[0];
	double csy = cs.y[0];
	double csz = cs.z[0];
	
	for(size_t j = 0; j < a.size2(); ++j) {
		p =  csx*q.x + csy*q.y + csz*q.z;
		double sign_sin = (p<0) ? -1.0 : 1.0;
		ap = abs(p);
		p = ap - long(ap/M_PI_twice)*M_PI_twice - M_PI;		// wrap p to -pi..pi		
		cp = ( p<M_PI_half )  ? p + M_PI_half : p - M_PI_3half;
		
		//pre-fetch next data
		csx = cs.x[j+1];
		csy = cs.y[j+1];
		csz = cs.z[j+1];	

		a(iframe,j) = complex<double>(sfs[j]*sign_sin*sine(p),sfs[j]*sine(cp));
	}	
}

void AllScatterDevice::scatter_frames(CartesianCoor3D& q) {
	
	for(size_t i = 0; i < myframes.size(); ++i)
	{
		scatter_frame(i,q);
	}
}

void AllScatterDevice::norm1() {
	// prepare scatter factors 
	vector<double> scatterfactors(a.size2());
	
	for(size_t i = 0; i < a.size1(); ++i)
	{
		for(size_t j = 1; j < a.size2(); ++j)
		{
			a(i,0) += a(i,j);
		}
	}

}

vector<complex<double> > AllScatterDevice::correlate_frames() {
	// each nodes has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)

	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());

	vector<complex<double> > local_A;
	local_A.resize(maxCSsize,complex<double>(0,0));
	
	for(size_t ci = 0; ci < myframes.size(); ++ci)
	{
		local_A[ci] = a(ci,0);
	}


	vector<double> local_Ar = flatten(local_A);
	vector<double> all_Ar; all_Ar.resize(2*maxCSsize*p_thisworldcomm->size());	

	boost::mpi::all_gather(*p_thisworldcomm,&local_Ar[0], 2*maxCSsize ,&all_Ar[0]);

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

	vector<double> ain_r = flatten(correlated_a);
	
	vector<double> aout_r; aout_r.resize(ain_r.size());
	boost::mpi::reduce(*p_thisworldcomm,&ain_r[0],ain_r.size(),&aout_r[0],std::plus<double>(),0);

	vector<complex<double> > aout = compress(aout_r);
	return aout;
}

vector<complex<double> > AllScatterDevice::conjmultiply_frames() {
	// each node has computed their assigned frames
	// the total scattering amplitudes reside in a(x,0)


	//negotiate maximum size for coordinatesets
	size_t CSsize = myframes.size();
	size_t maxCSsize;
	boost::mpi::all_reduce(*p_thisworldcomm,CSsize,maxCSsize,boost::mpi::maximum<size_t>());


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

void AllScatterDevice::superpose_spectrum(vector<complex<double> >& spectrum, vector<complex<double> >& fullspectrum) {
	for(size_t j = 0; j < spectrum.size(); ++j)
	{
		fullspectrum[j] += spectrum[j];
	}
}

void AllScatterDevice::execute(CartesianCoor3D& q) {
			
	string avm = Params::Inst()->scattering.average.method;
				
	VectorUnfold* p_vectorunfold = NULL;			
	if (avm=="bruteforce") {
		string avv = Params::Inst()->scattering.average.vectors;
		string avt = Params::Inst()->scattering.average.type;

		if (avt=="sphere") {
			SphereVectorUnfold* p_vu = new SphereVectorUnfold(q);
			p_vu->set_resolution(Params::Inst()->scattering.average.resolution);
			p_vu->set_vectors(Params::Inst()->scattering.average.vectors);
			p_vu->set_seed(0);
			p_vectorunfold = p_vu; 
		}
		if (avt=="cylinder") {
			CylinderVectorUnfold* p_vu = new CylinderVectorUnfold(q);
			p_vu->set_resolution(Params::Inst()->scattering.average.resolution);
			p_vu->set_vectors(Params::Inst()->scattering.average.vectors);
			p_vu->set_axis(Params::Inst()->scattering.average.axis);
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

	vector<CartesianCoor3D>& qvectors = p_vectorunfold->qvectors();

	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; spectrum.resize(p_sample->frames.size());

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	for(size_t qi = 0; qi < qvectors.size(); ++qi)
	{		
//		scatter_frames(qvectors[qi]);
//		norm1(); // put summed scattering amplitudes into first atom entry
		scatter_frames_norm1(qvectors[qi]); // put summed scattering amplitudes into first atom entry
	
		vector<complex<double> > thisspectrum;
		if (Params::Inst()->scattering.correlation.type=="time") {
			if (Params::Inst()->scattering.correlation.method=="direct") {
				thisspectrum = correlate_frames(); // if correlation, otherwise do a elementwise conj multiply here			
			} else {
				Err::Inst()->write("Correlation method not understood. Supported methods: direct");
				throw;
			}
		} else {
			thisspectrum = conjmultiply_frames();
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

vector<complex<double> >& AllScatterDevice::get_spectrum() {
	return m_spectrum;
}
