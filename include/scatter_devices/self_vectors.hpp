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
	boost::mpi::communicator m_scattercomm;
	boost::mpi::communicator m_fqtcomm;

	Sample* p_sample;
	
	std::vector< std::complex<double> >* p_asingle; 
	
	ScatterFactors scatterfactors;	
	
	std::vector<std::pair<size_t,CartesianCoor3D> > m_qvectorindexpairs;
	size_t m_current_qvector;
	bool m_writeflag;
	std::string m_fqt_filename;
	
	std::vector<std::complex<double> > m_spectrum;

    std::map<size_t,std::vector<CartesianCoor3D> > m_all_postalignmentvectors;
    
	std::vector<CartesianCoor3D> qvectors;
    
    std::vector<size_t> m_indexes;
    std::vector<std::vector<double> > m_x;
    std::vector<std::vector<double> > m_y;
    std::vector<std::vector<double> > m_z;
    
    
	void scatter(size_t ai, size_t mi);	
	
    void init(CartesianCoor3D& q);
	void correlate();
    void infinite_correlate();
    void average_correlate();    
    void norm();
    size_t get_numberofmoments();
    
    void multiply_alignmentfactors(size_t mi);
    void gather_sum();
		
public: 
	SelfVectorsScatterDevice(
			boost::mpi::communicator scatter_comm,
			boost::mpi::communicator fqt_comm,
			Sample& sample,
			std::vector<std::pair<size_t,CartesianCoor3D> > QIV,
			std::string fqt_filename
	);

	void compute();
	void next();
	void write();
	double progress();
};

#endif

//end of file
