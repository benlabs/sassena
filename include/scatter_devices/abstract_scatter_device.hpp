/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__ABSTRACT_SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__ABSTRACT_SCATTER_DEVICE_HPP_

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
#include <fftw3.h>

// other headers
#include "decomposition/assignment.hpp"
#include "math/coor3d.hpp"
#include "report/timer.hpp"
#include "services.hpp"
#include "scatter_devices/scatter_factors.hpp"



 template<typename Data>
 class concurrent_queue
 {
 private:
     std::queue<Data> the_queue;
     mutable boost::mutex the_mutex;
     boost::condition_variable the_condition_variable;
 public:

     void push(Data const& data)
     {
         boost::mutex::scoped_lock lock(the_mutex);
         the_queue.push(data);
         lock.unlock();

         // notify everyone who tries to push
         the_condition_variable.notify_all();         
     }

     bool empty() const
     {
         boost::mutex::scoped_lock lock(the_mutex);
         return the_queue.empty();
     }

     size_t size() const
     {
         return the_queue.size();
     }

     bool try_pop(Data& popped_value)
     {
         boost::mutex::scoped_lock lock(the_mutex);
         if(the_queue.empty())
         {
             return false;
         }

         popped_value=the_queue.front();
         the_queue.pop();
         // notify everyone who tries to push
         return true;
     }

     void wait_and_pop(Data& popped_value)
     {
         boost::mutex::scoped_lock lock(the_mutex);
         while (the_queue.empty())
         {
             the_condition_variable.wait(lock);
         }

         popped_value=the_queue.front();
         the_queue.pop();
         
     }

     void wait_for_empty() {
         while (!the_queue.empty()) boost::this_thread::sleep(boost::posix_time::milliseconds(25));
     }
 };
 
class IScatterDevice {
protected:
    virtual void runner() = 0;

    virtual size_t status() = 0;
    virtual double progress() = 0;

public:
    Timer timer;

    virtual void run() = 0;
};

class AbstractScatterDevice : public IScatterDevice {
protected:
    coor_t* p_coordinates;
    
    boost::mpi::communicator allcomm_;
    boost::mpi::communicator partitioncomm_;
	Sample& sample_;

	std::vector<CartesianCoor3D> vectors_;
    size_t current_vector_;
    
    boost::shared_ptr<MonitorClient> p_monitor_;
    boost::shared_ptr<HDF5WriterClient> p_hdf5writer_;
    
    size_t NN,NF,NA;
    
    fftw_complex* atfinal_;
    DivAssignment assignment_;
    
    ScatterFactors scatterfactors;
        
    virtual void stage_data() = 0;
    virtual void compute_serial() = 0;

    void next();
    void write();
    
    void runner();
    
    virtual void print_pre_stage_info() {}
    virtual void print_post_stage_info() {}
    virtual void print_pre_runner_info() {}
    virtual void print_post_runner_info() {}
    
    // use by threaded version
    virtual void compute_threaded() = 0;
    virtual void start_workers() = 0;
    virtual void stop_workers() = 0;
        
    
    size_t status();
    double progress();
    
public:

    AbstractScatterDevice(
        boost::mpi::communicator allcomm,
        boost::mpi::communicator partitioncomm,
        Sample& sample,
        std::vector<CartesianCoor3D> vectors,
        size_t NAF,
        boost::asio::ip::tcp::endpoint fileservice_endpoint,
		boost::asio::ip::tcp::endpoint monitorservice_endpoint
        );    
    ~AbstractScatterDevice();
    
    void run();
};

#endif

//end of file
