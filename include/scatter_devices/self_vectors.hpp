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
#include "math/coor3d.hpp"
#include "sample.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

#include "scatter_devices/scatter_device.hpp"

class SelfVectorsScatterDevice : public ScatterDevice {
private:
	boost::mpi::communicator* p_thisworldcomm;

	Sample* p_sample;
	
	std::vector< std::complex<double> >* p_asingle; 
	
	ScatterFactors scatterfactors;	
	
    std::vector<CartesianCoor3D> qvectors;
	
	std::vector<std::complex<double> > m_spectrum;

    std::map<size_t,std::vector<CartesianCoor3D> > m_all_postalignmentvectors;
    
    
    std::vector<size_t> m_indexes;
    std::vector<std::vector<double> > m_x;
    std::vector<std::vector<double> > m_y;
    std::vector<std::vector<double> > m_z;
    
    
	void scatter(size_t ai, size_t mi);	
	
    void init(CartesianCoor3D& q);
	void correlate();
    void infinite_correlate();
    void norm();
    size_t get_numberofmoments();
    
    void multiply_alignmentfactors(size_t mi);
    void gather_sum();
		
public: 
	SelfVectorsScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);

	void execute(CartesianCoor3D q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
	
};

#endif

//end of file
