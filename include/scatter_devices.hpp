/*
 *  scatterdevices.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTERDEVICES_HPP_
#define SCATTERDEVICES_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpi.hpp>

// other headers
#include "atomselection.hpp"
#include "coor3d.hpp"
#include "coordinate_sets.hpp"
#include "frame.hpp"
#include "frames.hpp"
#include "particle_trajectory.hpp"
#include "sample.hpp"
#include "scatter_factors.hpp"
#include "timer.hpp"

class ScatterDevice {

	Timer timer;	

public: 
	
	virtual void execute(CartesianCoor3D& q) = 0;
	virtual std::vector<std::complex<double> >& get_spectrum() = 0;
	
};

class SelfScatterDevice : public ScatterDevice {

	boost::mpi::communicator* p_thisworldcomm;

	Sample* p_sample;
	
	boost::numeric::ublas::matrix<std::complex<double> > a; // rows = particles, columns = time coordinate
	
	std::vector<ParticleTrajectory> particle_trajectories;
	
	ScatterFactors scatterfactors;	
	
	std::vector<std::complex<double> > m_spectrum;

	void scatter_particle(size_t iparticle, CartesianCoor3D& q);	
	void scatter_particles(CartesianCoor3D& q);	
	
	// special routines for no-time correlated calculation:
	void scatter_particle_instant(size_t iparticle, CartesianCoor3D& q);	
	void scatter_particles_instant(CartesianCoor3D& q);	
	
	
	void assemble_spectrum();
	void crosssum_particles();
	void correlate_particles();
	void conjmultiply_particles();
	
	void superpose_spectrum(std::vector<std::complex<double> >& fullspectrum);	

	
public: 
	SelfScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);

	void execute(CartesianCoor3D& q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
	
};


class AllScatterDevice : public ScatterDevice {

	boost::mpi::communicator* p_thisworldcomm;

	Sample* p_sample;

	boost::numeric::ublas::matrix<std::complex<double> > a; // rows = time coordinate, columns = particles
	
	CoordinateSets coordinate_sets;
	
	ScatterFactors scatterfactors;
		
	std::vector<size_t> myframes;
		
	std::vector<std::complex<double> > m_spectrum;
		
	void scatter_frame(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void norm1();
	
	void scatter_frame_norm1(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames_norm1(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude

	std::vector<std::complex<double> > correlate_frames();
	std::vector<std::complex<double> > conjmultiply_frames();

	
	void superpose_spectrum(std::vector<std::complex<double> >& spectrum, std::vector<std::complex<double> >& fullspectrum);	
	
public: 
	AllScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);
	
	void execute(CartesianCoor3D& q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
};


#endif

//end of file