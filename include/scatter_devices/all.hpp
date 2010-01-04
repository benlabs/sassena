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

#ifndef SCATTER_DEVICES__ALL_HPP_
#define SCATTER_DEVICES__ALL_HPP_

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
#include "sample.hpp"

#include "math/coor3d.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

#include "scatter_devices/scatter_device.hpp"

class AllScatterDevice : public ScatterDevice {
protected:
	boost::mpi::communicator* p_thisworldcomm;
	Sample* p_sample;

    // first = q, second = frames
	std::vector< std::vector< std::complex<double> > >* p_a; 
	std::vector< std::complex<double> >* p_asingle; 
	
	ScatterFactors scatterfactors;
	std::vector<size_t> myframes;
		
	std::vector<std::complex<double> > m_spectrum;		

	// have to be implemented by concrete classes:
    virtual void init(CartesianCoor3D& q) = 0;
    virtual void norm() = 0;
    virtual size_t get_numberofmoments() = 0 ;	
	virtual void scatter(size_t moffset,size_t mcount) = 0 ;


    void multiply_alignmentfactors(CartesianCoor3D q);
    
    void exchange();
    void correlate();
    void gather_sum();
    
    void conjmultiply();
    void sum();
    void gather_cat();
	
public: 
	AllScatterDevice(boost::mpi::communicator& thisworld, Sample& sample);
	~AllScatterDevice();
	
	void execute(CartesianCoor3D q); 
	std::vector<std::complex<double> >& get_spectrum(); // returns F(q,tau)
};


#endif

//end of file
