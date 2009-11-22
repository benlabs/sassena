/*
 *  scatter_devices/self_vectors.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__SELF_VECTORS_HPP_
#define SCATTER_DEVICES__SELF_VECTORS_HPP_

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
#include "coor3d.hpp"
#include "particle_trajectory.hpp"
#include "sample.hpp"
#include "scatter_factors.hpp"
#include "scatter_devices.hpp"
#include "timer.hpp"

#include "scatter_devices/scatter_device.hpp"

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

#endif

//end of file
