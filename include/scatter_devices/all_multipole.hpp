/*
 *  scatter_devices/all_multipole.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__ALL_MULTIPOLE_HPP_
#define SCATTER_DEVICES__ALL_MULTIPOLE_HPP_

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
#include "sample.hpp"
#include "particle_trajectory.hpp"
#include "scatter_factors.hpp"
#include "timer.hpp"

#include "scatter_devices/scatter_device.hpp"

// define a parent class for both spherical and clyindrial 
class AllMSScatterDevice : public ScatterDevice {

	boost::mpi::communicator* p_thisworldcomm;

	Sample* p_sample;

	boost::numeric::ublas::matrix<std::complex<double> > a; // rows = time coordinate, columns = particles

	ScatterFactors scatterfactors;

	std::vector<size_t> myframes;

	std::vector<std::complex<double> > m_spectrum;

	void scatter_frame(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void norm1();
	
	void scatter_frame_norm1(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames_norm1(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude

	std::vector<std::complex<double> > correlate_frames_fftw(long mpindex);
	std::vector<std::complex<double> > correlate_frames(long mpindex);
	void conjmultiply_frames();
    void multiply_alignmentfactors(CartesianCoor3D q);


	std::vector<std::complex<double> > gather_frames();

	void superpose_spectrum(std::vector<std::complex<double> >& spectrum, std::vector<std::complex<double> >& fullspectrum);	

public: 
	AllMSScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);

	void execute(CartesianCoor3D q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
};


// define a parent class for both spherical and clyindrial 
class AllMCScatterDevice : public ScatterDevice {

	boost::mpi::communicator* p_thisworldcomm;

	Sample* p_sample;

	boost::numeric::ublas::matrix<std::complex<double> > a; // rows = time coordinate, columns = multipole moments

	ScatterFactors scatterfactors;

	std::vector<size_t> myframes;

	std::vector<std::complex<double> > m_spectrum;

	void scatter_frame(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void norm1();

	void scatter_frame_norm1(size_t iframe,CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude
	void scatter_frames_norm1(CartesianCoor3D& q); // a(x,0) contains the total scattering amplitude

	std::vector<std::complex<double> > correlate_frames_fftw(long mpindex);
	std::vector<std::complex<double> > correlate_frames(long mpindex);
	void conjmultiply_frames();
    void multiply_alignmentfactors(CartesianCoor3D q);

	std::vector<std::complex<double> > gather_frames();

	void superpose_spectrum(std::vector<std::complex<double> >& spectrum, std::vector<std::complex<double> >& fullspectrum);	

public: 
	AllMCScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);

	void execute(CartesianCoor3D q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
};


#endif

// end of file

