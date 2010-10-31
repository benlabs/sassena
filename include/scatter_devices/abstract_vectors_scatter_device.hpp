/*
 *  scatter_devices/scatter_devices.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__ABSTRACT_VECTORS_SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__ABSTRACT_VECTORS_SCATTER_DEVICE_HPP_

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
#include "math/coor3d.hpp"
#include "report/timer.hpp"
#include "scatter_devices/abstract_scatter_device.hpp"

class AbstractVectorsScatterDevice : public AbstractScatterDevice {
protected:
    size_t NM;
    
	std::vector<CartesianCoor3D> subvector_index_;
    size_t current_subvector_;
        
    double progress();
    void init_subvectors(CartesianCoor3D& q);
    
    void print_pre_stage_info();
    void print_post_stage_info();
    void print_pre_runner_info();
    void print_post_runner_info();
        
public:
    AbstractVectorsScatterDevice(
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
