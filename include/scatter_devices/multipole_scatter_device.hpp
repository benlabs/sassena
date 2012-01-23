/** \file 
This file contains a class which implements the scattering calculation for multipole moment based orientally averaged scattering ( which is all type scattering ).

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/


#ifndef SCATTER_DEVICES__ABSTRACT_MULTIPOLE_SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__ABSTRACT_MULTIPOLE_SCATTER_DEVICE_HPP_

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
Implements all type scattering using the multipole method for spherical orientational averaging
*/
class MPSphereScatterDevice : public AbstractScatterDevice {
protected:
    
    // from abstract vectors scatter device
    size_t NM;
    long LMAX;
    
    CartesianCoor3D qvector_;
	std::vector<std::pair<long,long> > multipole_index_;
    size_t current_moment_;
    
    std::vector<std::complex<double> > j_l;
        
    double progress();
    void init_moments(CartesianCoor3D& q);
    
    void print_pre_stage_info();
    void print_post_stage_info();
    void print_pre_runner_info();
    void print_post_runner_info();
    
	bool ram_check();
	
    // from all vectors scatter device
    
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
    
    ~MPSphereScatterDevice();
    fftw_plan fftw_planF_;
    fftw_plan fftw_planB_;
    
        
public:
    MPSphereScatterDevice(
        boost::mpi::communicator allcomm,
        boost::mpi::communicator partitioncomm,
        Sample& sample,
        std::vector<CartesianCoor3D> vectors,
        size_t NAF,
        boost::asio::ip::tcp::endpoint fileservice_endpoint,
		boost::asio::ip::tcp::endpoint monitorservice_endpoint
        );
};

/** 
Implements all type scattering using the multipole method for cylindrical orientational averaging
*/
class MPCylinderScatterDevice : public AbstractScatterDevice {
protected:
    
    // from abstract vectors scatter device
    size_t NM;
        
    CartesianCoor3D qvector_;
	std::vector<std::pair<long,long> > multipole_index_;
    size_t current_moment_;
    
    std::vector<std::complex<double> > j_l;
        
    double progress();
    void init_moments(CartesianCoor3D& q);
    
    void print_pre_stage_info();
    void print_post_stage_info();
    void print_pre_runner_info();
    void print_post_runner_info();
    
	bool ram_check();
	
    // from all vectors scatter device
    
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
    
    ~MPCylinderScatterDevice();
    fftw_plan fftw_planF_;
    fftw_plan fftw_planB_;
    
        
public:
    MPCylinderScatterDevice(
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