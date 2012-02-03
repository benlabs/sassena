/** \file
This file contains a class which implements the scattering calculation for all scattering.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "stager/data_stager.hpp"

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

bool AllVectorsScatterDevice::ram_check() {
	// inherit ram requirements for parent class
	bool state = AbstractVectorsScatterDevice::ram_check();
	
    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
	size_t NMAXF = zeronode_assignment.max();
	size_t NNPP = partitioncomm_.size();
	size_t NTHREADS = Params::Inst()->limits.computation.threads;
    
	size_t memscale = Params::Inst()->limits.computation.memory.scale;
	// direct calculation buffer
	size_t bytesize_signal_buffer=0;
	size_t bytesize_exchange_buffer=0;
	size_t bytesize_alignpad_buffer=0;

	if (NNPP==1) {
		bytesize_signal_buffer=NF*NTHREADS*sizeof(fftw_complex);
		bytesize_exchange_buffer=0; // no exchange
		bytesize_alignpad_buffer=2*NF*sizeof(fftw_complex);
	} else {
		bytesize_signal_buffer=NMAXF*NNPP*sizeof(fftw_complex);
		bytesize_exchange_buffer=NMAXF*NNPP*sizeof(fftw_complex);
		bytesize_alignpad_buffer=2*NF*sizeof(fftw_complex);
	}

    if (bytesize_signal_buffer>memscale*Params::Inst()->limits.computation.memory.signal_buffer) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("limits.computation.memory.signal_buffer too small");
			Err::Inst()->write(string("limits.computation.memory.signal_buffer=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.signal_buffer));
			Err::Inst()->write(string("limits.computation.memory.scale=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.scale));	
            Err::Inst()->write(string("requested: ")+boost::lexical_cast<string>(bytesize_signal_buffer));            
        }
		state=false;
    }
    if (bytesize_exchange_buffer>memscale*Params::Inst()->limits.computation.memory.exchange_buffer) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("limits.computation.memory.exchange_buffer too small");
			Err::Inst()->write(string("limits.computation.memory.exchange_buffer=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.exchange_buffer));
			Err::Inst()->write(string("limits.computation.memory.scale=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.scale));	
            Err::Inst()->write(string("requested: ")+boost::lexical_cast<string>(bytesize_exchange_buffer));            
        }
		state=false;
    }
    if (bytesize_alignpad_buffer>memscale*Params::Inst()->limits.computation.memory.alignpad_buffer) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("limits.computation.memory.alignpad_buffer too small");
			Err::Inst()->write(string("limits.computation.memory.alignpad_buffer=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.alignpad_buffer));
			Err::Inst()->write(string("limits.computation.memory.scale=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.scale));				
            Err::Inst()->write(string("requested: ")+boost::lexical_cast<string>(bytesize_alignpad_buffer));            
        }
		state=false;
    }

	if (allcomm_.rank()==0) {
	Info::Inst()->write(string("required limits.computation.memory.signal_buffer=")+
		boost::lexical_cast<string>(bytesize_signal_buffer));
	Info::Inst()->write(string("required limits.computation.memory.exchange_buffer=")+
		boost::lexical_cast<string>(bytesize_exchange_buffer));
	Info::Inst()->write(string("required limits.computation.memory.alignpad_buffer=")+
		boost::lexical_cast<string>(bytesize_alignpad_buffer));
	}
	// final buffer for reduction:
	// NF*sizeof(fftw_complex)
	// however, at that moment the signal_buffer is deallocated
	// so we can ignore this requirement
	// since we're interested in the potential peak consumption..

	return state;
}

void AllVectorsScatterDevice::stage_data() {
    Timer& timer = timer_[boost::this_thread::get_id()];
    if (allcomm_.rank()==0) Info::Inst()->write(string("Forcing stager.mode=frames"));
    DataStagerByFrame data_stager(sample_,allcomm_,partitioncomm_,timer);
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

void AllVectorsScatterDevice::worker() {
      
    workerbarrier->wait();

    Timer& timer = timer_[boost::this_thread::get_id()];
     
    while (true) {
        size_t subvector_index;
        timer.start("sd:worker:wait");
        atscatter_.wait_and_pop(subvector_index);  
        timer.stop("sd:worker:wait");
        
        if (subvector_index<NM) {
            timer.start("sd:worker:scatter");
            scatter(subvector_index);
            timer.stop("sd:worker:scatter");
        }

        workerbarrier->wait();
    }
}

fftw_complex* AllVectorsScatterDevice::exchange() {
    
    size_t NNPP = partitioncomm_.size();
  
    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
  
    fftw_complex* atOUT = (fftw_complex*) fftw_malloc(NMAXF*NNPP*sizeof(fftw_complex)); 

    double* pIN = (double*) at_;
    double* pOUT = (double*) atOUT;

    boost::mpi::all_to_all(partitioncomm_,pIN,2*NMAXF,pOUT);

    return atOUT;
}

fftw_complex* AllVectorsScatterDevice::alignpad(fftw_complex* at) {

    size_t NNPP = partitioncomm_.size();

    DivAssignment zeronode_assignment(NNPP,0,NF);
    size_t NMAXF = zeronode_assignment.max();

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
    complex<double> a = smath::reduce<double>(at,NF) * (1.0/NF);
    afinal_ += a;
    a2final_ += a*conj(a);
    smath::add_elements(atfinal_,at,NF);
}

void AllVectorsScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    Timer& timer = timer_[boost::this_thread::get_id()];

    timer.start("sd:c:init");
    init_subvectors(q);    
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
    timer.stop("sd:c:init");

    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
    size_t NNPP = partitioncomm_.size();
            
    current_subvector_=0;
	memset(atfinal_,0,NF*sizeof(fftw_complex));
    afinal_ = 0;
    a2final_ = 0;
    
    timer.start("sd:c:block");
    // special case: 1 core, no exchange required
    if (NNPP==1) {
        size_t NTHREADS = worker_threads.size();
        
        at_ = (fftw_complex*) fftw_malloc(NF*NTHREADS*sizeof(fftw_complex));
        memset(at_,0,NF*NTHREADS*sizeof(fftw_complex));

        for(size_t i = 0; i < NM; i+=NTHREADS) {
            timer.start("sd:c:b:scatter");
            scatterblock(i,std::min(NTHREADS,NM-i));
            timer.stop("sd:c:b:scatter");
            
            timer.start("sd:c:b:dspstore");
            for(size_t j = 0; j < std::min(NTHREADS,NM-i); ++j)
            {
                size_t offset = ((j+i)%NTHREADS)*NF;
                fftw_complex* pat = &(at_[offset]);
                fftw_complex* nat = alignpad(pat);
                dsp(nat);
                store(nat);         
                fftw_free(nat); nat=NULL;
            }
            timer.stop("sd:c:b:dspstore");
            
            current_subvector_+=NTHREADS;
            p_monitor_->update(allcomm_.rank(),progress());
        }
        fftw_free(at_);        
        
    } else {
        at_ = (fftw_complex*) fftw_malloc(NMAXF*NNPP*sizeof(fftw_complex));
        memset(at_,0,NMAXF*NNPP*sizeof(fftw_complex));

        for(size_t i = 0; i < NM; i+=NNPP) {
            timer.start("sd:c:b:scatter");
            scatterblock(i,std::min(NNPP,NM-i));
            timer.stop("sd:c:b:scatter");
            
            timer.start("sd:c:b:exchange");
            fftw_complex* at = exchange();
            timer.stop("sd:c:b:exchange");
                        
            timer.start("sd:c:b:dspstore");
            if (partitioncomm_.rank()<std::min(NNPP,NM-i)) {
                fftw_complex* nat = alignpad(at);
                fftw_free(at); at=NULL;
                dsp(nat);
                store(nat);         
                fftw_free(nat); nat=NULL;  
            }
            if (at!=NULL) fftw_free(at); 
            timer.stop("sd:c:b:dspstore");
            
            current_subvector_+=std::min(NNPP,NM-i);
            timer.start("sd:c:b:progress");
            p_monitor_->update(allcomm_.rank(),progress());    	
            timer.stop("sd:c:b:progress");
        }
        fftw_free(at_);
    }
    timer.stop("sd:c:block");

    timer.start("sd:c:wait");
    partitioncomm_.barrier();
    timer.stop("sd:c:wait");
    
    timer.start("sd:c:reduce");
    if (NNPP>1) {

        double* p_atfinal = (double*) &(atfinal_[0][0]);
        double* p_atlocal=NULL;
        if (partitioncomm_.rank()==0) {
            fftw_complex* atlocal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
            p_atlocal = (double*) &(atlocal_[0][0]);
            memset(atlocal_,0,NF*sizeof(fftw_complex));        
        }

        boost::mpi::reduce(partitioncomm_,p_atfinal,2*NF,p_atlocal,std::plus<double>(),0);
        double* p_afinal = (double*) &(afinal_);
        std::complex<double> alocal(0);
        double* p_alocal = (double*) &(alocal);
        boost::mpi::reduce(partitioncomm_,p_afinal,2,p_alocal,std::plus<double>(),0);
        double* p_a2final = (double*) &(a2final_);
        std::complex<double> a2local(0);
        double* p_a2local = (double*) &(a2local);
        boost::mpi::reduce(partitioncomm_,p_a2final,2,p_a2local,std::plus<double>(),0);

        // swap pointers on rank 0 
        if (partitioncomm_.rank()==0) {
            fftw_free(atfinal_); 
            atfinal_ = (fftw_complex*) p_atlocal;
            afinal_ = alocal;
			a2final_ = a2local;
        }
    }
    timer.stop("sd:c:reduce");
    
    double factor = 1.0/subvector_index_.size();
    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
        afinal_ *= factor;        
        a2final_ *= factor;                
    }
}

void AllVectorsScatterDevice::scatterblock(size_t index,size_t count) {
    size_t SCATTERBLOCK = worker_threads.size();
    
//    cerr << "loop w/ index=" << index << " , count=" << count << endl;
//    sleep(3);
    for(size_t i = 0; i < count; ++i)
    {
        atscatter_.push(index+i);
        
        if (((i+1)%SCATTERBLOCK)==0) {
            workerbarrier->wait();
        }
    }
    
    if ((count%SCATTERBLOCK)!=0) {
        size_t rest = SCATTERBLOCK - (count%SCATTERBLOCK);
        for(size_t i = 0; i < rest; ++i)
        {
            atscatter_.push(NM); // NM will be ignored, but triggers a barrier            
        }
        workerbarrier->wait();
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
   size_t NMAXF = zeronode_assignment.max();
   size_t NNPP = partitioncomm_.size();
   size_t offset = 0;
   if (NNPP==1) {
       size_t NTHREADS = worker_threads.size();
       offset = (this_subvector%NTHREADS)*NF;
   } else {
       offset = (this_subvector%NNPP)*NMAXF;       
   }

   size_t NMYF = DivAssignment(partitioncomm_.size(),partitioncomm_.rank(),NF).size();
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
           Ar += esf*cos(p);
           Ai += esf*sin(p);
	   }
       p_a[fi][0] = Ar;
       p_a[fi][1] = Ai;
	}
}

// end of file
