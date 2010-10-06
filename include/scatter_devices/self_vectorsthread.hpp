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

#ifndef SCATTER_DEVICES__SELF__VECTORSTHREAD_HPP_
#define SCATTER_DEVICES__SELF__VECTORSTHREAD_HPP_

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
#include <boost/asio.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpi.hpp>
#include <boost/thread.hpp>

// other headers
#include "sample.hpp"
#include "services.hpp"
#include "math/coor3d.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

#include "scatter_devices/scatter_device.hpp"

class SelfVectorsThreadScatterDevice : public ScatterDevice {
protected:
    
    //////////////////////////////
    // data
    //////////////////////////////
    
	boost::mpi::communicator m_scattercomm;
	boost::mpi::communicator m_fqtcomm;
	Sample* p_sample;

	std::vector<std::pair<size_t,CartesianCoor3D> > m_qvectorindexpairs;
	size_t m_current_qvector;
    
    HDF5WriterClient* p_hdf5writer;
    MonitorClient* p_monitor;    
    
    size_t NN,NF,NA,NM;

    // first = q, second = frames
    concurrent_queue< std::pair<size_t,size_t> > at0;    
    concurrent_queue< std::vector< std::complex<double> >* > at1;
    std::vector< std::complex<double> > atfinal;

	std::vector<std::complex<double> > m_spectrum;		

    std::queue<boost::thread*> worker_threads;

    mutable boost::mutex at2_mutex;
    coor_t* p_coordinates;
    
	ScatterFactors scatterfactors;
	std::vector<size_t> m_assignment;
    
    std::vector<CartesianCoor3D> qvectors;
    
    //////////////////////////////
    // methods
    //////////////////////////////
    
    // computation
    std::vector<std::complex<double> >* scatter(size_t qindex,size_t aindex);
    void init(CartesianCoor3D& q);
	
    
    // organization
    void stage_data();

        
	void compute(bool marshal);
	void next();
	void write();
	
    void runner();
    
    void start_workers();
    void stop_workers();
    void worker1(bool loop);
    void worker2(bool loop);
    
    volatile size_t worker2_counter;
    mutable boost::mutex worker2_mutex;    
    volatile bool worker2_done;
    boost::condition_variable worker2_notifier;
    
public:
    SelfVectorsThreadScatterDevice(
			boost::mpi::communicator scatter_comm,
			boost::mpi::communicator fqt_comm,
			Sample& sample,
			std::vector<std::pair<size_t,CartesianCoor3D> > QVI,
			std::vector<size_t> assignment,
			boost::asio::ip::tcp::endpoint fileservice_endpoint,
			boost::asio::ip::tcp::endpoint monitorservice_endpoint			
	);
    virtual ~SelfVectorsThreadScatterDevice();

	size_t status();
    double progress();
    
    void run();
    
};

#endif