/** \file 
This file contains a class which generates the scattering device based on user input provides through the CONTROL module. The constructing routine only return a pointer to the interface.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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

/** 
Factory class which constructs the proper scattering device based on parameter settings.
*/
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
