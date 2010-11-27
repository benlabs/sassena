/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
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
#include <boost/thread/barrier.hpp>

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
    Timer timer;
    DataStagerByFrame data_stager(sample_,allcomm_,partitioncomm_,assignment_,timer);
    p_coordinates = data_stager.stage();
}

AllVectorsScatterDevice::~AllVectorsScatterDevice() {
    if (p_coordinates!=NULL) {
        free(p_coordinates);
        p_coordinates = NULL;
    }
    fftw_destroy_plan(fftw_planF_);
    fftw_destroy_plan(fftw_planB_);
    if (atfinal_!=NULL) {
        fftw_free(atfinal_);    
        atfinal_=NULL;
    }
}

void AllVectorsScatterDevice::start_workers() {

    //size_t worker3_threads = 1; // fix this at the moment;
    size_t worker1_threads = Params::Inst()->limits.computation.threads.scatter;

    // don't allow more worker1 threads than NM are available!
    if (worker1_threads>NM) worker1_threads = NM;
    
    if (allcomm_.rank()==0) {
        if (partitioncomm_.size()<worker1_threads) {
            Warn::Inst()->write("Partition smaller than number of threads.");
            Warn::Inst()->write("Can not utilize more threads than partition size.");
        }
    }
    
    for(size_t i = 0; i < worker1_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&AllVectorsScatterDevice::worker_scatter,this) ));	
    }
    
    scatterbarrier = new boost::barrier(worker1_threads+1);
}

void AllVectorsScatterDevice::stop_workers() {

    // implement hangups , i.e. sentinal values

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
    delete scatterbarrier;
}

void AllVectorsScatterDevice::worker_scatter() {
            
    while (true) {
        size_t subvector_index;
        atscatter_.wait_and_pop(subvector_index);  
//        cerr << "popped " << subvector_index << endl;
//        sleep(3);
        if (subvector_index<NM) {
            scatter(subvector_index);
        }

        scatterbarrier->wait();
    }
}

fftw_complex* AllVectorsScatterDevice::exchange() {
    
    size_t NNPP = partitioncomm_.size();
  
    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.size();
  
    fftw_complex* atOUT = (fftw_complex*) fftw_malloc(NMAXF*NNPP*sizeof(fftw_complex)); 

    double* pIN = (double*) at_;
    double* pOUT = (double*) atOUT;

    boost::mpi::all_to_all(partitioncomm_,pIN,2*NMAXF,pOUT);

    return atOUT;
}

fftw_complex* AllVectorsScatterDevice::alignpad(fftw_complex* at) {

    size_t NNPP = partitioncomm_.size();

    DivAssignment zeronode_assignment(NNPP,0,NF);
    size_t NMAXF = zeronode_assignment.size();

    fftw_complex* atOUT = (fftw_complex*) fftw_malloc(2*NF*sizeof(fftw_complex)); 

    for(size_t i = 0; i < NNPP; ++i)
    {
        DivAssignment node_assignment(NNPP,i,NF); 
        size_t offset = node_assignment.offset();
        size_t len = node_assignment.size();
        double* pIN = (double*) &(at[i*NMAXF]);
        double* pOUT = (double*) &(atOUT[offset]);
        memcpy(pOUT,pIN,len*sizeof(fftw_complex));      
    }
    memset(atOUT[NF],0,NF*sizeof(fftw_complex));

    return atOUT;
}

void AllVectorsScatterDevice::dsp(fftw_complex* at) {
    
	// correlate or sum up
    if (Params::Inst()->scattering.dsp.type=="autocorrelate") {
        if (Params::Inst()->scattering.dsp.method=="direct") {
            smath::auto_correlate_direct(at,NF);	        
        } else if (Params::Inst()->scattering.dsp.method=="fftw") {
            smath::auto_correlate_fftw(at,fftw_planF_,fftw_planB_,NF);
        } else {
            Err::Inst()->write("Correlation method not understood");
            Err::Inst()->write("scattering.dsp.method == direct, fftw");        
            throw;
        }
	} else if (Params::Inst()->scattering.dsp.type=="square")  {
        smath::square_elements(at,NF);
	} else if (!(Params::Inst()->scattering.dsp.type=="plain")) {
        Err::Inst()->write(string("DSP type not understood: ")+Params::Inst()->scattering.dsp.type);
        Err::Inst()->write("scattering.dsp.type == autocorrelate, square, plain");        
        throw;	    
    }
}

void AllVectorsScatterDevice::store(fftw_complex* at) {
    smath::add_elements(atfinal_,at,NF);
}

void AllVectorsScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    init_subvectors(q);
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop

    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.size();
    size_t NNPP = partitioncomm_.size();
            
    current_subvector_=0;
	memset(atfinal_,0,NF*sizeof(fftw_complex));

    size_t NMBLOCK = NNPP;
    
    // special case: 1 core, no exchange required
    if (NMBLOCK==1) {
        size_t NTHREADS = worker_threads.size();
        at_ = (fftw_complex*) fftw_malloc(NF*NTHREADS*sizeof(fftw_complex));
        memset(at_,0,NF*NTHREADS*sizeof(fftw_complex));

        for(size_t i = 0; i < NM; i+=NTHREADS) {
            scatterblock(i,std::min(NTHREADS,NM-i));
            for(size_t j = 0; j < std::min(NTHREADS,NM-i); ++j)
            {
                size_t offset = (j%NTHREADS)*NF;
                fftw_complex* pat = &(at_[offset]);
                fftw_complex* nat = alignpad(pat);
                dsp(nat);
                store(nat);         
                fftw_free(nat); nat=NULL;
            }
        }
                       
    } else {
        at_ = (fftw_complex*) fftw_malloc(NMAXF*NMBLOCK*sizeof(fftw_complex));
        memset(at_,0,NMAXF*NMBLOCK*sizeof(fftw_complex));

        for(size_t i = 0; i < NM; i+=NMBLOCK) {
            scatterblock(i,std::min(NMBLOCK,NM-i));

            fftw_complex* at = exchange();
            if (partitioncomm_.rank()<std::min(NMBLOCK,NM-i)) {
                fftw_complex* nat = alignpad(at);
                fftw_free(at); at=NULL;
                dsp(nat);
                store(nat);         
                fftw_free(nat); nat=NULL;  
            }
            if (at!=NULL) fftw_free(at); 
        }
        fftw_free(at_);
    }

    p_monitor_->update(allcomm_.rank(),progress());    	

    partitioncomm_.barrier();
    
    
    if (NN>1) {

        double* p_atfinal = (double*) &(atfinal_[0][0]);
        double* p_atlocal=NULL;
        if (partitioncomm_.rank()==0) {
            fftw_complex* atlocal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
            p_atlocal = (double*) &(atlocal_[0][0]);;
            memset(atlocal_,0,NF*sizeof(fftw_complex));        
        }

        boost::mpi::reduce(partitioncomm_,p_atfinal,2*NF,p_atlocal,std::plus<double>(),0);

        // swap pointers on rank 0 
        if (partitioncomm_.rank()==0) {
            fftw_free(atfinal_); 
            atfinal_ = (fftw_complex*) p_atlocal;
        }
    }
    
    double factor = 1.0/subvector_index_.size();
    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
    }
}

void AllVectorsScatterDevice::scatterblock(size_t index,size_t count) {
    size_t SCATTERBLOCK = Params::Inst()->limits.computation.threads.scatter;
    
//    cerr << "loop w/ index=" << index << " , count=" << count << endl;
//    sleep(3);
    for(size_t i = 0; i < count; ++i)
    {
        atscatter_.push(index+i);
        
        if (((i+1)%SCATTERBLOCK)==0) {
            scatterbarrier->wait();
        }
    }
    
    if ((count%SCATTERBLOCK)!=0) {
        size_t rest = SCATTERBLOCK - (count%SCATTERBLOCK);
        for(size_t i = 0; i < rest; ++i)
        {
            atscatter_.push(NM); // NM will be ignored, but triggers a barrier            
        }
        scatterbarrier->wait();
    }
}


void AllVectorsScatterDevice::scatter(size_t this_subvector) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;
//   cerr << "this_subvector " << this_subvector << endl;
//   sleep(3);
   std::vector<double>& sfs = scatterfactors.get_all();

   DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
   size_t NMAXF = zeronode_assignment.size();
   size_t NNPP = partitioncomm_.size();
   size_t NMBLOCK = NNPP;
   size_t offset = 0;
   if (NMBLOCK==1) {
       size_t NTHREADS = worker_threads.size();
       offset = (this_subvector%NTHREADS)*NF;
   } else {
       offset = (this_subvector%NMBLOCK)*NMAXF;       
   }

   size_t NMYF = assignment_.size();
   fftw_complex* p_a = &(at_[offset]);
   
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
}

// end of file
