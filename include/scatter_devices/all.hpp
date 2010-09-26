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
	boost::mpi::communicator m_scattercomm;
	boost::mpi::communicator m_fqtcomm;
	Sample* p_sample;


	std::vector<std::pair<size_t,CartesianCoor3D> > m_qvectorindexpairs;
	size_t m_current_qvector;
	bool m_writeflag;
	std::string m_fqt_filename;


    // first = q, second = frames
	std::vector< std::vector< std::complex<double> > >* p_a; 
	std::vector< std::complex<double> >* p_asingle; 
	
	ScatterFactors scatterfactors;
	std::vector<size_t> myframes;
    std::vector<CoordinateSet*> csets;
		
	std::vector<std::complex<double> > m_spectrum;		
	std::vector<CartesianCoor3D> qvectors;

	// have to be implemented by concrete classes:
    virtual void init(CartesianCoor3D& q) = 0;
    virtual void norm() = 0;
    virtual size_t get_numberofmoments() = 0 ;	
	virtual void scatter(size_t moffset,size_t mcount) = 0 ;


    void multiply_alignmentfactors(CartesianCoor3D q);
    
    void exchange();
    void correlate();
    void infinite_correlate();
    void average_correlate();
    void gather_sum();
    
    void conjmultiply();
    void sum();
    void gather_cat();
	
public: 
    AllScatterDevice(
			boost::mpi::communicator scatter_comm,
			boost::mpi::communicator fqt_comm,
			Sample& sample,
			std::vector<std::pair<size_t,CartesianCoor3D> > QVI,
			std::string fqt_filename
	);
    virtual ~AllScatterDevice();
	
	void compute();
	void next();
	void write();
	size_t status();
    double progress();
	
};


#endif

//end of file
