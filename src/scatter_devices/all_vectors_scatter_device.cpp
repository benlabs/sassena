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
#include "scatter_devices/all_vectors_scatter_device.hpp"

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

AllVectorsScatterDevice::AllVectorsScatterDevice(
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
//    fftw_free(wspace);
}

void AllVectorsScatterDevice::stage_data() {
    DataStagerByFrame data_stager(sample_,allcomm_,partitioncomm_,assignment_);
    p_coordinates = data_stager.stage();
}

AllVectorsScatterDevice::~AllVectorsScatterDevice() {
    if (p_coordinates!=NULL) {
        free(p_coordinates);
        p_coordinates = NULL;
    }
    fftw_destroy_plan(fftw_planF_);
    fftw_destroy_plan(fftw_planB_);
    if (at3_!=NULL) {
        fftw_free(at3_); 
        at3_=NULL;
    }
    if (atfinal_!=NULL) {
        fftw_free(atfinal_);    
        atfinal_=NULL;
    }
}

void AllVectorsScatterDevice::start_workers() {

    //size_t worker3_threads = 1; // fix this at the moment;
    size_t worker3_threads = Params::Inst()->limits.computation.worker3_threads;    
    size_t worker2_threads = 1;
    size_t worker1_threads = Params::Inst()->limits.computation.worker1_threads;

    // don't allow more worker1 threads than NM are available!
    if (worker1_threads>NM) worker1_threads = NM;
    
    
    for(size_t i = 0; i < worker1_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&AllVectorsScatterDevice::worker1,this) ));	
    }

    worker_threads.push( new boost::thread(boost::bind(&AllVectorsScatterDevice::worker2,this) ));	

    for(size_t i = 0; i < worker3_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&AllVectorsScatterDevice::worker3,this) ));	
    }
}

void AllVectorsScatterDevice::stop_workers() {

    // implement hangups , i.e. sentinal values

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
}


void AllVectorsScatterDevice::worker1_task(size_t subvector_index) {
    std::pair<size_t,fftw_complex* > p_a;
    p_a.first = subvector_index;
    p_a.second = scatter(subvector_index);
    at1_.push( p_a );
}

void AllVectorsScatterDevice::worker1() {
    
    size_t maxbuffer_bytesize = Params::Inst()->limits.memory.at1_buffer;
    size_t timeout = Params::Inst()->limits.times.at1_buffer;

    size_t single_entry_bytesize = assignment_.size()*sizeof(fftw_complex);
    size_t max_entries = maxbuffer_bytesize/single_entry_bytesize;
    
    while (true) {
        
        // retrieve from queue
        size_t subvector_index;

        // QoS here
        while (at1_.size()>max_entries) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(timeout));
        }
        at0_.wait_and_pop(subvector_index);  
        
        worker1_task(subvector_index);
        
    }
}

void AllVectorsScatterDevice::worker2_task(const size_t subvector_index,fftw_complex* p_a) {
    // compute global assignment patterns:
    
    vector<DivAssignment> assignments;
    for (size_t i=0;i<NN;i++) {
        assignments.push_back(DivAssignment(NN,i,NF));
    }
	
    size_t target_node = (subvector_index%NN);
    size_t rank = partitioncomm_.rank();
        
    if (rank==target_node) {
        // allocate an array of size 2*NF, to allow direct fftw operations
        fftw_complex* p_at_allframes = (fftw_complex*) fftw_malloc(2*NF*sizeof(fftw_complex));
        for(size_t i = 0; i < NN; ++i)
        {
            double* pp_at_allframes = (double*) &(p_at_allframes[assignments[i].offset()]);
            if (i==rank) {
                memcpy(pp_at_allframes,&(p_a[0][0]),assignments[i].size()*sizeof(fftw_complex));
            } else {
                partitioncomm_.recv(i,0,pp_at_allframes,2*assignments[i].size()); 
            }
        }   
        // zero pad the data before pushing
        memset(&p_at_allframes[NF],0,NF*sizeof(fftw_complex));
        
        at2_.push(p_at_allframes);  
    } else {
        double* p_at_frames = (double*) &(p_a[0][0]);
        partitioncomm_.send(target_node,0,p_at_frames,2*assignment_.size());
    }
    worker2_counter++;
    fftw_free(p_a); p_a=NULL;
}


void AllVectorsScatterDevice::worker2() {
    
    std::map<size_t,fftw_complex* > local_queue;
    
    while (true) {
        // retrieve from queue
        std::pair<size_t,fftw_complex*> at_local_pair(0,NULL);

        // first check local queue
        size_t key = worker2_counter;
        if (local_queue.find(key)==local_queue.end()) {
            at1_.wait_and_pop(at_local_pair);   
        } else {
            at_local_pair.first = key;
            at_local_pair.second = local_queue[key];
            local_queue.erase(key);
        }

        if (at_local_pair.first!=key) {
            local_queue[at_local_pair.first]=at_local_pair.second;
            continue;
        }
        worker2_task(at_local_pair.first,at_local_pair.second);

        if (worker2_counter==NM) {
            at2_.push(NULL);
        }                
        current_subvector_++;
        p_monitor_->update(allcomm_.rank(),progress());    						            
    }
    
}

void AllVectorsScatterDevice::worker3_task(fftw_complex* p_a) {
    
	// correlate or sum up
    if (Params::Inst()->scattering.dsp.type=="autocorrelate") {
        if (Params::Inst()->scattering.dsp.method=="direct") {
            smath::auto_correlate_direct(p_a,NF);	        
        } else if (Params::Inst()->scattering.dsp.method=="fftw") {
            smath::auto_correlate_fftw(p_a,fftw_planF_,fftw_planB_,NF);
        } else {
            Err::Inst()->write("Correlation method not understood");
            Err::Inst()->write("scattering.dsp.method == direct, fftw");        
            throw;
        }
	} else if (Params::Inst()->scattering.dsp.type=="square")  {
        smath::square_elements(p_a,NF);
	}
    boost::mutex::scoped_lock at3l(at3_mutex);
    smath::add_elements(at3_,p_a,NF);
    at3l.unlock();
    fftw_free(p_a); p_a=NULL;
}

void AllVectorsScatterDevice::worker3() {
    
    while ( true ) {
        fftw_complex* p_at_local=NULL;

        at2_.wait_and_pop(p_at_local);
        if (p_at_local==NULL) { // ats as sentinal
            worker3_done = true;
            boost::mutex::scoped_lock w3l(worker3_mutex);
            w3l.unlock();
            worker3_notifier.notify_all();
                       
            p_monitor_->update(allcomm_.rank(),progress());    						            
            continue;
        }
        
        worker3_task(p_at_local);
	}	    
	
}

void AllVectorsScatterDevice::compute_serial() {
	CartesianCoor3D q=vectors_[current_vector_];

    init_subvectors(q);

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
            
    current_subvector_=0;
    
    at3_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
    memset(at3_,0,NF*sizeof(fftw_complex));


    for(size_t i = 0; i < NM; ++i)
    {
        worker1_task(i);
        current_subvector_++;
        
        std::pair<size_t,fftw_complex* > at_local_pair(0,NULL);
        at1_.wait_and_pop(at_local_pair);
        worker2_task(at_local_pair.first,at_local_pair.second);
        fftw_complex* at_local = NULL;        
        at2_.try_pop(at_local);
        if (at_local!=NULL) worker3_task(at_local);
        p_monitor_->update(allcomm_.rank(),progress());
    }

    if (partitioncomm_.rank()==0) {
        atfinal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
        memset(atfinal_,0,NF*sizeof(fftw_complex));        
    }

    double factor = 1.0/subvector_index_.size();
    if (NN>1) {

        double* p_at3 = (double*) &(at3_[0][0]);
        double* p_atlocal=NULL;
        if (partitioncomm_.rank()==0) {
            p_atlocal = (double*) &(atfinal_[0][0]);;
        }

        boost::mpi::reduce(partitioncomm_,p_at3,2*NF,p_atlocal,std::plus<double>(),0);
    }
    else {
        memcpy(atfinal_,at3_,NF*sizeof(fftw_complex));
    }
    
    fftw_free(at3_); at3_=NULL;

    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
    }
}


void AllVectorsScatterDevice::compute_threaded() {
    CartesianCoor3D q=vectors_[current_vector_];

    init_subvectors(q);

	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
            
    current_subvector_=0;
	
    worker2_counter = 0;
    worker3_done = false;
    boost::mutex::scoped_lock w3l(worker3_mutex);

    at3_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
    memset(at3_,0,NF*sizeof(fftw_complex));

    for(size_t i = 0; i < NM; ++i) at0_.push(i);
    while (!worker3_done)  worker3_notifier.wait(w3l);

    if (partitioncomm_.rank()==0) {
        atfinal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
        memset(atfinal_,0,NF*sizeof(fftw_complex));        
    }

    double factor = 1.0/subvector_index_.size();
    if (NN>1) {

        double* p_at3 = (double*) &(at3_[0][0]);
        double* p_atlocal=NULL;
        if (partitioncomm_.rank()==0) {
            p_atlocal = (double*) &(atfinal_[0][0]);;
        }

        boost::mpi::reduce(partitioncomm_,p_at3,2*NF,p_atlocal,std::plus<double>(),0);
    }
    else {
        memcpy(atfinal_,at3_,NF*sizeof(fftw_complex));
    }
    
    fftw_free(at3_); at3_=NULL;

    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
    }
}

fftw_complex* AllVectorsScatterDevice::scatter(size_t this_subvector) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;
   
   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = assignment_.size();
   fftw_complex* p_a = (fftw_complex*) fftw_malloc(NMYF*sizeof(fftw_complex));
   
   CartesianCoor3D q = subvector_index_[this_subvector];
   double qx = q.x;
   double qy = q.y;
   double qz = q.z;
   coor_t x,y,z;
   double esf,p;

   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       coor_t* p_data = &(p_coordinates[fi*NA*3]);
               
       double Ar =0; 
       double Ai = 0;
    
	   for(size_t j = 0; j < NA; ++j) {   		

	       esf = sfs[j];
	       x = p_data[3*j];
 	       y = p_data[3*j+1];
	       z = p_data[3*j+2];
       
           p =  x*qx + y*qy + z*qz;
           
           Ar += sfs[j]*cos(p);
           Ai += sfs[j]*sin(p);
	   }
       p_a[fi][0] = Ar;
       p_a[fi][1] = Ai;
	}
	
    return p_a;
}

