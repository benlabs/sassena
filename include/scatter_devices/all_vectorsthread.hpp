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

#ifndef SCATTER_DEVICES__ALL__VECTORSTHREAD_HPP_
#define SCATTER_DEVICES__ALL__VECTORSTHREAD_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>
#include <queue>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpi.hpp>
#include <boost/thread.hpp>

// other headers
#include "sample.hpp"

#include "math/coor3d.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

#include "scatter_devices/scatter_device.hpp"

class AllVectorsThreadScatterDevice : public ScatterDevice {
protected:
	boost::mpi::communicator* p_thisworldcomm;
	Sample* p_sample;

    // first = q, second = frames
    std::queue< std::vector< std::complex<double> > > at1;
    boost::mutex at1_mutex;

    std::queue< std::vector< std::complex<double> > > at2;
    boost::mutex at2_mutex;

    std::vector< std::complex<double> > at3;
	
	ScatterFactors scatterfactors;
	std::vector<size_t> myframes;
		
	std::vector<std::complex<double> > m_spectrum;		

    void multiply_alignmentfactors(CartesianCoor3D q);
    
    void correlate(std::vector< std::complex<double> >& data);
    void infinite_correlate(std::vector< std::complex<double> >& data);
    void conjmultiply(std::vector< std::complex<double> >& data);
    
    std::vector<CartesianCoor3D> qvectors;
    
	// have to be implemented by concrete classes:
    void init(CartesianCoor3D& q);
    size_t get_numberofmoments();	
	void scatter(size_t moffset);
    
    void worker1(size_t mi, CartesianCoor3D q);
    void worker2();
    void worker3();
    
    bool worker1_done_flag;
    bool worker2_done_flag;
    bool worker3_done_flag;
	
public: 
	AllVectorsThreadScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);
	~AllVectorsThreadScatterDevice();
	
	void execute(CartesianCoor3D q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
};


#endif

//end of file
