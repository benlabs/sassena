/** \file 
This file contains an refined version of the abstract scatter device, used for performing vector based orientationally averaged scattering calculations.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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

/** 
Implements control flow for vector based scattering calculations 
*/
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

	bool ram_check();
        
public:
    AbstractVectorsScatterDevice(
        boost::mpi::communicator allcomm,
        boost::mpi::communicator partitioncomm,
        Sample& sample,
        std::vector<CartesianCoor3D> vectors,
        size_t NAF,
        boost::asio::ip::tcp::endpoint fileservice_endpoint,
		boost::asio::ip::tcp::endpoint monitorservice_endpoint
        );
};

#endif

//end of file
