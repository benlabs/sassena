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

#ifndef SCATTER_DEVICES__ALL__VECTORSTHREAD_HPP_
#define SCATTER_DEVICES__ALL__VECTORSTHREAD_HPP_

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

// other headers
#include "sample.hpp"

#include "math/coor3d.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

#include "scatter_devices/scatter_device.hpp"


 template<typename Data>
 class concurrent_queue
 {
 private:
     std::queue<Data> the_queue;
     mutable boost::mutex the_mutex;
     boost::condition_variable the_condition_variable;
     size_t max_size;
 public:
     void set_maxsize(size_t size) {
         max_size = size;
     }
     size_t get_maxsize(size_t size) {
         return max_size;
     }
     
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
 
class AllVectorsThreadScatterDevice : public ScatterDevice {
protected:
	boost::mpi::communicator m_scattercomm;
	boost::mpi::communicator m_fqtcomm;
	Sample* p_sample;

	std::vector<std::pair<size_t,CartesianCoor3D> > m_qvectorindexpairs;
	size_t m_current_qvector;
	bool m_writeflag;
	std::string m_fqt_filename;
    boost::asio::ip::tcp::endpoint m_fileserver_endpoint;

    // first = q, second = frames
    concurrent_queue< std::pair<size_t,std::vector< std::complex<double> >* > > at1;
    concurrent_queue< std::vector< std::complex<double> >* > at2;
    std::vector< std::complex<double> > at3;
	
    mutable boost::mutex scatter_mutex;
	
	ScatterFactors scatterfactors;
	std::vector<size_t> myframes;
    std::vector<CoordinateSet*> csets;
		
	std::vector<std::complex<double> > m_spectrum;		
    
    void multiply_alignmentfactors(CartesianCoor3D q);
    
    void correlate(std::vector< std::complex<double> >& data);
    void infinite_correlate(std::vector< std::complex<double> >& data);
    void conjmultiply(std::vector< std::complex<double> >& data);
    
    std::vector<CartesianCoor3D> qvectors;
    
	// have to be implemented by concrete classes:
    void init(CartesianCoor3D& q);
    size_t get_numberofmoments();	
	void scatter(size_t moffset);
    
    void worker1(size_t mi, CartesianCoor3D q);
    void worker2();
    void worker3();
    
    size_t worker1counter;
    size_t worker2counter;
    size_t worker3counter;
    size_t worker3maxcounter;
	
public: 
    AllVectorsThreadScatterDevice(
			boost::mpi::communicator scatter_comm,
			boost::mpi::communicator fqt_comm,
			Sample& sample,
			std::vector<std::pair<size_t,CartesianCoor3D> > QVI,
			boost::asio::ip::tcp::endpoint filemutex_server
	);
    virtual ~AllVectorsThreadScatterDevice();

	void compute();
	void next();
	void write();
	size_t status();
    double progress();
};


#endif

//end of file
