/*
 *  scatterdevices.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/self_vectors_scatter_device.hpp"

// standard header
#include <algorithm>
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "math/coor3d.hpp"
#include "math/smath.hpp"
#include "decomposition/decompose.hpp"
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "scatter_devices/data_stager.hpp"

using namespace std;

SelfVectorsScatterDevice::SelfVectorsScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<CartesianCoor3D> vectors,
    size_t NAF,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    AbstractVectorsScatterDevice(
        allcomm,
        partitioncomm,
        sample,
        vectors,
        NAF,
        fileservice_endpoint,
        monitorservice_endpoint
    )
{
    fftw_complex* wspace= NULL;
    fftw_planF_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_planB_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);    
}

void SelfVectorsScatterDevice::stage_data() {
    DataStagerByAtom data_stager(sample_,allcomm_,partitioncomm_,assignment_);
    p_coordinates = data_stager.stage();
}

SelfVectorsScatterDevice::~SelfVectorsScatterDevice() {
    if (p_coordinates!=NULL) {
        free(p_coordinates);
        p_coordinates=NULL;
    }
    fftw_destroy_plan(fftw_planF_);
    fftw_destroy_plan(fftw_planB_);
    if (at2_!=NULL) {
        fftw_free(at2_);
        at2_=NULL;
    }
    if (atfinal_!=NULL) {
        fftw_free(atfinal_);
        atfinal_=NULL;
    }
}

void SelfVectorsScatterDevice::start_workers() {

    size_t worker2_threads = 1;
    size_t worker1_threads = Params::Inst()->limits.computation.worker1_threads;
    
    // don't allow more worker1 threads than NM are available!
    if (worker1_threads>NM) worker1_threads = NM;
    
    
    for(size_t i = 0; i < worker1_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&SelfVectorsScatterDevice::worker1,this) ));	
    }

    worker_threads.push( new boost::thread(boost::bind(&SelfVectorsScatterDevice::worker2,this) ));	
}

void SelfVectorsScatterDevice::stop_workers() {

    // implement hangups , i.e. sentinal values

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
}

void SelfVectorsScatterDevice::worker1_task(size_t this_subvector,size_t this_atom) {
    fftw_complex* p_at_local = NULL;
    p_at_local=scatter(this_subvector,this_atom);
    
    // apply dsp
    if (Params::Inst()->scattering.dsp.type=="autocorrelate") {
        if (Params::Inst()->scattering.dsp.method=="direct") {
            smath::auto_correlate_direct(p_at_local,NF);	        
        } else if (Params::Inst()->scattering.dsp.method=="fftw") {          
            smath::auto_correlate_fftw(p_at_local,fftw_planF_,fftw_planB_,NF);
        } else {
            Err::Inst()->write("Correlation method not understood");
            Err::Inst()->write("scattering.dsp.method == direct, fftw");        
            throw;
        }
	} else if (Params::Inst()->scattering.dsp.type=="square") {
        smath::square_elements(p_at_local,NF);
	}
    at1_.push(p_at_local);
}
                
void SelfVectorsScatterDevice::worker1() {
    
    size_t maxbuffer_bytesize = Params::Inst()->limits.memory.at1_buffer;
    size_t timeout = Params::Inst()->limits.times.at1_buffer;

    size_t single_entry_bytesize = NF*sizeof(fftw_complex);
    size_t max_entries = maxbuffer_bytesize/single_entry_bytesize;
    
    while (true) {
        
        // retrieve from queue
        std::pair<size_t,size_t> qmoment_atom;

        // QoS here
        while (at1_.size()>max_entries) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(timeout));
        }
        at0_.wait_and_pop(qmoment_atom);

        worker1_task(qmoment_atom.first,qmoment_atom.second);

    }
}

void SelfVectorsScatterDevice::worker2_task(fftw_complex* p_a) {

    worker2_counter++;
    smath::add_elements(at2_,p_a,NF);
    p_monitor_->update(allcomm_.rank(),progress());    						            

    fftw_free(p_a); p_a=NULL;
}

void SelfVectorsScatterDevice::worker2() {
    while ( true ) {
        fftw_complex* p_at_local=NULL;
        at1_.wait_and_pop(p_at_local);
        if (p_at_local==NULL) { // ats as sentinal
            worker2_done = true;
            boost::mutex::scoped_lock w2l(worker2_mutex);
            w2l.unlock();
            worker2_notifier.notify_all();
                   
           p_monitor_->update(allcomm_.rank(),progress());    						            
           continue;
        }
        worker2_task(p_at_local);
        worker2_counter++;
    }
}

void SelfVectorsScatterDevice::compute_serial() {
	CartesianCoor3D q=vectors_[current_vector_];
    init_subvectors(q);

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop

    current_subvector_=0;

    at2_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
    memset(at2_,0,NF*sizeof(fftw_complex));

    for(size_t i = 0; i < NM; ++i)
    {
        for(size_t n = 0; n < assignment_.size(); ++n)
        {
            worker1_task(i,n);
            fftw_complex* p_a=NULL;        
            at1_.wait_and_pop(p_a);
            worker2_task(p_a); 
        }
        current_subvector_++;
        p_monitor_->update(allcomm_.rank(),progress());                                  
    }

    if (partitioncomm_.rank()==0) {
        atfinal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
        memset(atfinal_,0,NF*sizeof(fftw_complex));        
    }
    
    double factor = 1.0/subvector_index_.size();    
    if (NN>1) {
        double* p_at2 = (double*) &(at2_[0][0]);

        double* p_atlocal = NULL;
        if (partitioncomm_.rank()==0) {
            p_atlocal = (double*) &(atfinal_[0][0]);
        }

        boost::mpi::reduce(partitioncomm_,p_at2,2*NF,p_atlocal,std::plus<double>(),0);
    }
    else {
        memcpy(atfinal_,at2_,NF*sizeof(fftw_complex));
    }
    
    fftw_free(at2_); at2_=NULL;

    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
    }
}

void SelfVectorsScatterDevice::compute_threaded() {
    CartesianCoor3D q=vectors_[current_vector_];

    init_subvectors(q);

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
            
    current_subvector_=0;

    worker2_counter = 0;
    worker2_done = false;
    boost::mutex::scoped_lock w2l(worker2_mutex);

    at2_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
    memset(at2_,0,NF*sizeof(fftw_complex));


    for(size_t i = 0; i < NM; ++i) {
        for(size_t n = 0; n < assignment_.size(); ++n)
        {
            at0_.push(make_pair(i,n));
        }
    }
    while (!worker2_done)  worker2_notifier.wait(w2l);

    if (partitioncomm_.rank()==0) {
        atfinal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
        memset(atfinal_,0,NF*sizeof(fftw_complex));        
    }
    
    double factor = 1.0/subvector_index_.size();    
    if (NN>1) {
        double* p_at2 = (double*) &(at2_[0][0]);

        double* p_atlocal = NULL;
        if (partitioncomm_.rank()==0) {
            p_atlocal = (double*) &(atfinal_[0][0]);
        }

        boost::mpi::reduce(partitioncomm_,p_at2,2*NF,p_atlocal,std::plus<double>(),0);
    }
    else {
        memcpy(atfinal_,at2_,NF*sizeof(fftw_complex));
    }
    
    fftw_free(at2_); at2_=NULL;

    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
    }
}

fftw_complex* SelfVectorsScatterDevice::scatter(size_t mi,size_t ai) {

	// this is broken <-- revise this!!!
	double s = scatterfactors.get(assignment_[ai]);
	
	// double allocate (2*NF), to allow direct application of autocorrelation. 
    fftw_complex* p_at_local = (fftw_complex*) fftw_malloc(2*NF*sizeof(fftw_complex));
    
    double qx = subvector_index_[mi].x;
    double qy = subvector_index_[mi].y;
    double qz = subvector_index_[mi].z;
    
    coor_t* p_data = &(p_coordinates[ai*NF*3]);
	
    	
	for(size_t j = 0; j < NF; ++j)
	{
		coor_t x1 = p_data[j*3];
		coor_t y1 = p_data[j*3 + 1];
		coor_t z1 = p_data[j*3 + 2];

		double p1 = x1*qx+y1*qy+z1*qz;
        double sp1 = sin(p1);
        double cp1 = cos(p1);
        p_at_local[j][0] = s*cp1;
        p_at_local[j][1] = s*sp1;
	}
	
    memset(&p_at_local[NF],0,NF*sizeof(fftw_complex));
  
    return p_at_local;
}
