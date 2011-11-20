/** \file 
This file contains a class which implements the scattering calculation for all scattering.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include <boost/thread/barrier.hpp>
#include <fftw3.h>

// other headers
#include "sample.hpp"
#include "math/coor3d.hpp"
#include "report/timer.hpp"

#include "scatter_devices/abstract_vectors_scatter_device.hpp"

/** 
Implements all type scattering using vectors for orientational averaging
*/
class AllVectorsScatterDevice : public AbstractVectorsScatterDevice {
protected:

    // first = q, second = frames
    fftw_complex* at_;
        
	// data, outer loop by frame, inner by atoms, XYZ entries
    coor_t* p_coordinates;
    
	void scatter(size_t this_subvector);
    
    void stage_data();
   
    concurrent_queue< size_t > atscatter_;    
    void worker();        
	void compute();	

    void scatterblock(size_t index,size_t count);
    void store(fftw_complex* at);
    void dsp(fftw_complex* at);
    fftw_complex* alignpad(fftw_complex* at);
    fftw_complex* exchange();
    
    ~AllVectorsScatterDevice();
    fftw_plan fftw_planF_;
    fftw_plan fftw_planB_;
public: 
    AllVectorsScatterDevice(
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
