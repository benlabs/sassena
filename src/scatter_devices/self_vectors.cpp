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
#include "control.hpp"
#include "sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"
#include "vector_unfold.hpp"


using namespace std;

SelfScatterDevice::SelfScatterDevice(boost::mpi::communicator& thisworld, Sample& sample){

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;
	
	size_t nn = thisworld.size(); // Number of Nodes
	size_t na = sample.atoms.selections[target].size(); // Number of Atoms
	size_t nf = sample.coordinate_sets.size();

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
	
	timer.start("sd:init:scatter:assign");	
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
	timer.stop("sd:init:scatter:assign");	

	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

	timer.start("sd:init:scatter:data");	
	for(size_t r = 0; r < nn; ++r)
	{

		vector<size_t> assigned_frames = edecomp.indexes_for(r);
		
		for(size_t i = 0; i < assigned_frames.size(); ++i)
		{
			double  *p_xc,*p_yc,*p_zc;	
			
			CoordinateSet* p_cs =NULL;
			if (rank==r) {
                p_cs = &( p_sample->coordinate_sets.load(assigned_frames[i]) );
                                
				p_xc = &(p_cs->c1[0]);
				p_yc = &(p_cs->c2[0]);
				p_zc = &(p_cs->c3[0]);
							
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
				p_xc = &(p_cs->c1[minatoms]);
				p_yc = &(p_cs->c2[minatoms]);
				p_zc = &(p_cs->c3[minatoms]);
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
				
		}
	}
	timer.stop("sd:init:scatter:data");	
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
	a.resize(particle_trajectories.size(),sample.coordinate_sets.size());

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
	
}

void SelfScatterDevice::scatter_particles(CartesianCoor3D& q) {
	for(size_t i = 0; i < a.size1(); ++i)
	{
		timer.start("sd:fs:f");	    
		scatter_particle(i,q);
		timer.stop("sd:fs:f");		
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
		timer.stop("sd:fs:f");		    
		scatter_particle_instant(i,q);
		timer.stop("sd:fs:f");			
	}
}

void SelfScatterDevice::correlate_particles() {

    vector<complex<double> > acopy(a.size2());
    	
	for(size_t i = 0; i < a.size1(); ++i)
	{
	    
	    for(size_t j = 0; j < a.size2(); ++j) // make a working copy for the current particle i
	    {
            acopy[j] = a(i,j);
	    }
	    
	    // overwrite a with correlated values
	    for(size_t j = 0; j < a.size2(); ++j) // this iterates tau
	    {
            a(i,j)=0;

		    size_t last_starting_frame = a.size2()-j;
		    for(size_t k = 0; k < last_starting_frame; ++k) // this iterates the starting frame
		    {
			    a(i,j) += ( acopy[k] )*conj( acopy[j+k] ) ;
		    }

    	    a(i,j) *= (1.0/last_starting_frame); 
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
    
   string avm = Params::Inst()->scattering.average.orientation.method;
   
   VectorUnfold* p_vectorunfold = NULL;			
   timer.start("sd:vector:unfold"); 
   if (avm=="vectors") {
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
   
   vector<CartesianCoor3D>& qvectors = p_vectorunfold->vectors();
   timer.stop("sd:vector:unfold"); 
	
	
	/// k, qvectors are prepared:
	vector<complex<double> > spectrum; 
	
	if (p_thisworldcomm->rank()==0)  spectrum.resize(a.size2()); // only rank0 node aggregrates results

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");	
	
	for(size_t qi = 0; qi < qvectors.size(); ++qi)
	{
		vector<complex<double> > thisspectrum;
		if (Params::Inst()->scattering.correlation.type=="time") {
		    timer.start("sd:fs");
			scatter_particles(qvectors[qi]);
		    timer.stop("sd:fs");
		    
			if (Params::Inst()->scattering.correlation.method=="direct") {
			    timer.start("sd:correlate");					
				correlate_particles(); // if correlation, otherwise do a elementwise conj multiply here	
				timer.stop("sd:correlate");					
				
			} else {
				Err::Inst()->write("Correlation method not understood. Supported methods: direct");
				throw;
			}
		} else {
		    timer.start("sd:fs");		    
			scatter_particles_instant(qvectors[qi]); // this version takes advantage of a*a(conj) = (real,0)
		    timer.stop("sd:fs");
		    			
			timer.start("sd:conjmul");			    			
			conjmultiply_particles();
			timer.stop("sd:conjmul");				
		}
		crosssum_particles(); // the spectrum is in a(0,x)		
		assemble_spectrum(); // now rank==0 has the full spectrum in a(0,x)


		if (p_thisworldcomm->rank()==0) {
			timer.start("sd:superpose");		    
		    superpose_spectrum(spectrum);
			timer.stop("sd:superpose");
	    }
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

// end of file
