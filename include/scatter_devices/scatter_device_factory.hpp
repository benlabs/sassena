/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__SCATTER_DEVICE_FACTORY_HPP_
#define SCATTER_DEVICES__SCATTER_DEVICE_FACTORY_HPP_

// common header
#include "common.hpp"

#include <boost/asio.hpp>
#include <boost/mpi.hpp>

// other headers
#include "sample.hpp"
#include "scatter_devices/abstract_scatter_device.hpp"

class ScatterDeviceFactory {
public: 
    static IScatterDevice* create(
    		boost::mpi::communicator& scatter_comm,
    		Sample& sample,
    		boost::asio::ip::tcp::endpoint fileservice_endpoint,
    		boost::asio::ip::tcp::endpoint monitorservice_endpoint,
    		std::vector<CartesianCoor3D>& qvectors);
    	
};


#endif

//end of file
