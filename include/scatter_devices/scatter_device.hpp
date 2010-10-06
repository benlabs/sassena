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

#ifndef SCATTER_DEVICES__SCATTER_DEVICE_HPP_
#define SCATTER_DEVICES__SCATTER_DEVICE_HPP_

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
 
class ScatterDevice {
    
    virtual void runner() = 0;

    virtual size_t status() = 0;
    virtual double progress() = 0;

public:
    Timer timer;

//    ScatterDevice (
//        boost::mpi::communicator all_comm,
//	    boost::mpi::communicator partition_comm,
//	    Sample& sample,
//	    vector<pair<size_t,CartesianCoor3D> > QIV,
//	    boost::asio::ip::tcp::endpoint fileservice_endpoint,
//        boost::asio::ip::tcp::endpoint monitorservice_endpoint
//    );
    
//    Timer& get_timer() { return timer; }

    virtual void run() = 0;
};

//class AVScatterDevice : public ScatterDevice {
//    void stage_data();
//    
//    
//public:
//    AVScatterDevice (
//        boost::mpi::communicator all_comm,
//	    boost::mpi::communicator partition_comm,
//	    Sample& sample,
//	    vector<pair<size_t,CartesianCoor3D> > QIV,
//	    boost::asio::ip::tcp::endpoint fileservice_endpoint,
//        boost::asio::ip::tcp::endpoint monitorservice_endpoint
//    );
//    
//    
//};
//
//

#endif

//end of file
