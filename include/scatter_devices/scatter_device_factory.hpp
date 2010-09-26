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

#ifndef SCATTER_DEVICES__SCATTER_DEVICE_FACTORY_HPP_
#define SCATTER_DEVICES__SCATTER_DEVICE_FACTORY_HPP_

// common header
#include "common.hpp"


#include <boost/mpi.hpp>

// other headers
#include "sample.hpp"
#include "scatter_devices/scatter_device.hpp"

class ScatterDeviceFactory {
public: 
    static ScatterDevice* create(
    		boost::mpi::communicator& scatter_comm,
    		Sample& sample,
    		std::vector<CartesianCoor3D>& qvectors);
    	
};


#endif

//end of file
