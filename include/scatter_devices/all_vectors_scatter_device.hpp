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

#ifndef SCATTER_DEVICES__ALL_VECTORS_SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__ALL_VECTORS_SCATTER_DEVICE_HPP_

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
#include <fftw3.h>

// other headers
#include "sample.hpp"
#include "math/coor3d.hpp"
#include "report/timer.hpp"

#include "scatter_devices/abstract_vectors_scatter_device.hpp"

class AllVectorsScatterDevice : public AbstractVectorsScatterDevice {
protected:

    // first = q, second = frames
    concurrent_queue< size_t > at0_;    
    concurrent_queue< std::pair<size_t,std::vector< std::complex<double> >* > > at1_;
    concurrent_queue< std::vector< std::complex<double> >* > at2_;
    mutable std::vector< std::complex<double> > at3_;
    mutable boost::mutex at3_mutex;
	std::queue<boost::thread*> worker_threads;
    
	// data, outer loop by frame, inner by atoms, XYZ entries
    coor_t* p_coordinates;
    
	std::vector< std::complex<double> >* scatter(size_t this_subvector);
    
    void stage_data();
        
    void worker1();
    void worker2();
    void worker3();

    void worker1_task(size_t this_subvector);
    void worker2_task(size_t this_subvector,std::vector<std::complex<double> >* p_a);
    void worker3_task(std::vector<std::complex<double> >* p_a);
    
    volatile size_t worker2_counter;
    mutable boost::mutex worker3_mutex;
    mutable boost::mutex worker2_mutex;    
    volatile bool worker2_done;
    volatile bool worker3_done;
    boost::condition_variable worker3_notifier;
    boost::condition_variable worker2_notifier;
    
	void compute_serial();
	void compute_threaded();	
	void start_workers();
    void stop_workers();


    ~AllVectorsScatterDevice();
    fftw_plan fftw_planF_;
    fftw_plan fftw_planB_;
public: 
    AllVectorsScatterDevice(
			boost::mpi::communicator allcomm,
			boost::mpi::communicator partitioncomm,
			Sample& sample,
			std::vector<CartesianCoor3D> vectors,
			std::vector<size_t> assignment,
			boost::asio::ip::tcp::endpoint fileservice_endpoint,
			boost::asio::ip::tcp::endpoint monitorservice_endpoint			

	);
};


#endif

//end of file
