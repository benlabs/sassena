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
    atfinal_ = NULL;
    
    fftw_complex* wspace= NULL;
    fftw_planF_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_planB_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);    
}

void SelfVectorsScatterDevice::stage_data() {
    Timer timer;
    DataStagerByAtom data_stager(sample_,allcomm_,partitioncomm_,assignment_,timer);
    p_coordinates = data_stager.stage();
}

SelfVectorsScatterDevice::~SelfVectorsScatterDevice() {
    if (p_coordinates!=NULL) {
        free(p_coordinates);
        p_coordinates=NULL;
    }
    fftw_destroy_plan(fftw_planF_);
    fftw_destroy_plan(fftw_planB_);
    if (atfinal_!=NULL) {
        fftw_free(atfinal_);
        atfinal_=NULL;
    }
}

void SelfVectorsScatterDevice::start_workers() {

    size_t worker1_threads = Params::Inst()->limits.computation.threads.scatter;
    
    // don't allow more worker1 threads than NM are available!
    if (worker1_threads>NM) worker1_threads = NM;
    
    
    for(size_t i = 0; i < worker1_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&SelfVectorsScatterDevice::worker_scatter,this) ));	
    }
    
    scatterbarrier = new boost::barrier(worker1_threads+1);
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

void SelfVectorsScatterDevice::dsp(fftw_complex* at) {
    
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
                
void SelfVectorsScatterDevice::store(fftw_complex* at) {
    smath::add_elements(atfinal_,at,NF);
}


void SelfVectorsScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    init_subvectors(q);
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
            
    current_subvector_=0;
	memset(atfinal_,0,NF*sizeof(fftw_complex));
    
    size_t NTHREADS = worker_threads.size();
    // double allocate for zero padding
    at_ = (fftw_complex*) fftw_malloc(2*NF*NTHREADS*sizeof(fftw_complex));
    memset(at_,0,2*NF*NTHREADS*sizeof(fftw_complex));

    for(size_t n = 0; n < assignment_.size(); ++n)
    {
        for(size_t i = 0; i < NM; i+=NTHREADS) {
            scatterblock(n,i,std::min(NTHREADS,NM-i));
            for(size_t j = 0; j < std::min(NTHREADS,NM-i); ++j)
            {
                size_t offset = ((j+i)%NTHREADS)*2*NF;
                fftw_complex* pat = &(at_[offset]);
                dsp(pat);
                store(pat);         
            }
        }
    }

    p_monitor_->update(allcomm_.rank(),progress());    	

    partitioncomm_.barrier();
    
//    timer.start("sd:reduce");
    if (NN>1) {
        double* p_atfinal = (double*) &(atfinal_[0][0]);

        double* p_atlocal = NULL;
        if (partitioncomm_.rank()==0) {
            fftw_complex* atlocal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));            
            p_atlocal = (double*) &(atlocal_[0][0]);
            memset(p_atlocal,0,NF*sizeof(fftw_complex));        
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
//    timer.stop("sd:reduce");
    
}


void SelfVectorsScatterDevice::scatterblock(size_t atomindex, size_t index,size_t count) {
    size_t SCATTERBLOCK = Params::Inst()->limits.computation.threads.scatter;
    
//    cerr << "loop w/ index=" << index << " , count=" << count << endl;
//    sleep(3);
    for(size_t i = 0; i < count; ++i)
    {
        atscatter_.push(make_pair(index+i,atomindex));
        
        if (((i+1)%SCATTERBLOCK)==0) {
            scatterbarrier->wait();
        }
    }
    
    if ((count%SCATTERBLOCK)!=0) {
        size_t rest = SCATTERBLOCK - (count%SCATTERBLOCK);
        for(size_t i = 0; i < rest; ++i)
        {
            atscatter_.push(make_pair(NM,atomindex)); // NM will be ignored, but triggers a barrier            
        }
        scatterbarrier->wait();
    }
}

void SelfVectorsScatterDevice::worker_scatter() {
            
    while (true) {
        std::pair<size_t,size_t> vectoratom;
        atscatter_.wait_and_pop(vectoratom);  
//        cerr << "popped " << subvector_index << endl;
//        sleep(3);
        if (vectoratom.first<NM) {
            scatter(vectoratom.first,vectoratom.second);
        }

        scatterbarrier->wait();
    }
}

fftw_complex* SelfVectorsScatterDevice::scatter(size_t mi,size_t ai) {

	// this is broken <-- revise this!!!
	double s = scatterfactors.get(assignment_[ai]);
	
	// double allocate (2*NF), to allow direct application of autocorrelation. 
	size_t NTHREADS = worker_threads.size();
    size_t offset = (mi%NTHREADS)*2*NF;
    
    fftw_complex* p_at_local = &(at_[offset]);
    
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

// end of file
