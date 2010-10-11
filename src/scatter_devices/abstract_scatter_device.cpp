/*
 *  scatter_device.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/abstract_scatter_device.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "scatter_devices/data_stager.hpp"

using namespace std;

AbstractScatterDevice::AbstractScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<std::pair<size_t,CartesianCoor3D> > vector_index,
    std::vector<size_t> assignment,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),
    sample_(sample),
    vector_index_(vector_index),
    current_vector_(0),
    assignment_(assignment)
{
    
    p_hdf5writer_ = boost::shared_ptr<HDF5WriterClient>(new HDF5WriterClient(fileservice_endpoint));
    p_monitor_ = boost::shared_ptr<MonitorClient>(new MonitorClient(monitorservice_endpoint));
    	
    NN = partitioncomm_.size();
	NF = sample_.coordinate_sets.size();
	std::string target = Params::Inst()->scattering.target;
	NA = sample_.atoms.selections[target].indexes.size(); // Number of Atoms
	
	sample_.coordinate_sets.set_selection(sample.atoms.selections[target]);
	
	scatterfactors.set_sample(sample_);
	scatterfactors.set_selection(sample_.atoms.selections[target]);
	scatterfactors.set_background(true);
	
    atfinal_.resize(NF);
}

void AbstractScatterDevice::run() {
    
    print_pre_stage_info();
    timer.start("sd:stage");
    stage_data();
    timer.stop("sd:stage");
    
    print_post_stage_info();
    

    print_pre_runner_info();
    timer.start("sd:runner");
    runner();
    timer.stop("sd:runner");
    print_post_runner_info();

    // notify the services to finalize
    p_monitor_->update(allcomm_.rank(),1.0);
    p_hdf5writer_->flush();        
}


void AbstractScatterDevice::runner() {
 
    bool threads_on = Params::Inst()->limits.computation.threads;
    
    if (threads_on) {
         start_workers();
         while(status()==0) {
         	compute_threaded();
     		write();		    
            next();
         }             
         stop_workers();
    } else {
        while(status()==0) {
            compute_serial();
    		write();		    
            next();
        }                     
    }
}


void AbstractScatterDevice::next() {
	if (current_vector_>=vector_index_.size()) return;
	current_vector_++;    
}

double AbstractScatterDevice::progress() {
    double scale = 1.0/vector_index_.size();
    double base =  current_vector_*scale;
    return base;    
}

size_t AbstractScatterDevice::status() {
    if (current_vector_==vector_index_.size()) return 1; else return 0;
}

void AbstractScatterDevice::write() {
    if (partitioncomm_.rank()==0) {
        size_t index = vector_index_[current_vector_].first; 
        p_hdf5writer_->write(index,atfinal_);
    }
}