/** \file
This file contains a class which implements the scattering calculation for self scattering.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "stager/data_stager.hpp"

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
    ),
    assignment_(partitioncomm_.size(),partitioncomm_.rank(),NAF)
{    
    atfinal_ = NULL;
	
    fftw_complex* wspace= NULL;
    fftw_planF_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_planB_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);    
}

void SelfVectorsScatterDevice::stage_data() {
    Timer& timer = timer_[boost::this_thread::get_id()];   
    if (allcomm_.rank()==0) Info::Inst()->write(string("Forcing stager.mode=atoms")); 
    DataStagerByAtom data_stager(sample_,allcomm_,partitioncomm_,timer);
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

bool SelfVectorsScatterDevice::ram_check() {
	// inherit ram requirements for parent class
	bool state = AbstractVectorsScatterDevice::ram_check();
	
	// total memory requirements during computation:
    // atfinal = NF
    // at_ = NTHREADS*2*NF
    
	size_t NTHREADS = Params::Inst()->limits.computation.threads;

	size_t bytesize_signal_buffer=0;
	size_t memscale = Params::Inst()->limits.computation.memory.scale;
    
	bytesize_signal_buffer=2*NF*NTHREADS*sizeof(fftw_complex);
	
	if (bytesize_signal_buffer>memscale*Params::Inst()->limits.computation.memory.signal_buffer) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("limits.computation.memory.signal_buffer too small");
			Err::Inst()->write(string("limits.computation.memory.signal_buffer=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.signal_buffer));
			Err::Inst()->write(string("limits.computation.memory.scale=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.scale));			
            Err::Inst()->write(string("requested: ")+boost::lexical_cast<string>(bytesize_signal_buffer));            
        }
		state=false;
    }
    
	// final buffer for reduction:
	// NF*sizeof(fftw_complex)
	// however, at that moment the signal_buffer is deallocated
	// so we can ignore this requirement
	// since we're interested in the potential peak consumption..
	if (allcomm_.rank()==0) {
	Info::Inst()->write(string("required limits.computation.memory.signal_buffer=")+
		boost::lexical_cast<string>(bytesize_signal_buffer));
	}
	return state;
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
    complex<double> a = smath::reduce<double>(at,NF) * (1.0/NF);
    afinal_ += a;
    a2final_ += a*conj(a);  
    smath::add_elements(atfinal_,at,NF);
}


void SelfVectorsScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    Timer& timer = timer_[boost::this_thread::get_id()];

    timer.start("sd:c:init");
    init_subvectors(q);
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
    timer.stop("sd:c:init");
            
    current_subvector_=0;
    current_atomindex_=0;
	memset(atfinal_,0,NF*sizeof(fftw_complex));
    afinal_=0; 
    a2final_=0;
    
    size_t NTHREADS = worker_threads.size();
	size_t NNPP = partitioncomm_.size();
    at_ = (fftw_complex*) fftw_malloc(2*NF*NTHREADS*sizeof(fftw_complex));
    memset(at_,0,2*NF*NTHREADS*sizeof(fftw_complex));

    timer.start("sd:c:block");
    for(size_t n = 0; n < assignment_.size(); ++n)
    {
        current_subvector_=0;
        for(size_t i = 0; i < NM; i+=NTHREADS) {
            timer.start("sd:c:b:scatter");
            scatterblock(n,i,std::min(NTHREADS,NM-i));
            timer.stop("sd:c:b:scatter");

            timer.start("sd:c:b:dspstore");
            for(size_t j = 0; j < std::min(NTHREADS,NM-i); ++j)
            {
                size_t offset = ((j+i)%NTHREADS)*2*NF;
                fftw_complex* pat = &(at_[offset]);
                dsp(pat);
                store(pat);         
            }
            timer.stop("sd:c:b:dspstore");
            
            current_subvector_+=std::min(NTHREADS,NM-i);

            timer.start("sd:c:b:progress");
            p_monitor_->update(allcomm_.rank(),progress());
            timer.stop("sd:c:b:progress");

        }
        current_atomindex_+=1;
    }
    timer.stop("sd:c:block");
    
    fftw_free(at_);
    
    timer.start("sd:c:wait");
    partitioncomm_.barrier();
    timer.stop("sd:c:wait");
    
    timer.start("sd:c:reduce");
    if (NNPP>1) {
        double* p_atfinal = (double*) &(atfinal_[0][0]);

        double* p_atlocal = NULL;
        if (partitioncomm_.rank()==0) {
            fftw_complex* atlocal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));            
            p_atlocal = (double*) &(atlocal_[0][0]);
            memset(p_atlocal,0,NF*sizeof(fftw_complex));        
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


void SelfVectorsScatterDevice::scatterblock(size_t atomindex, size_t index,size_t count) {
    size_t SCATTERBLOCK = worker_threads.size();
    
//    cerr << "loop w/ index=" << index << " , count=" << count << endl;
//    sleep(3);
    for(size_t i = 0; i < count; ++i)
    {
        atscatter_.push(std::make_pair(index+i,atomindex));
        
        if (((i+1)%SCATTERBLOCK)==0) {
            workerbarrier->wait();
        }
    }
    
    if ((count%SCATTERBLOCK)!=0) {
        size_t rest = SCATTERBLOCK - (count%SCATTERBLOCK);
        for(size_t i = 0; i < rest; ++i)
        {
            atscatter_.push(std::make_pair(NM,atomindex)); // NM will be ignored, but triggers a barrier            
        }
        workerbarrier->wait();
    }
}

void SelfVectorsScatterDevice::worker() {

    workerbarrier->wait();

    Timer& timer = timer_[boost::this_thread::get_id()];
                
    while (true) {
        std::pair<size_t,size_t> vectoratom;
        timer.start("sd:worker:wait");        
        atscatter_.wait_and_pop(vectoratom);  
        timer.stop("sd:worker:wait");
        
        if (vectoratom.first<NM) {
            timer.start("sd:worker:scatter");
            scatter(vectoratom.first,vectoratom.second);
            timer.stop("sd:worker:scatter");
        }

        workerbarrier->wait();
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


double SelfVectorsScatterDevice::progress() {
    double scale1 = 1.0/vectors_.size();
    double scale2 = 1.0/assignment_.size();
    double scale3 = 1.0/subvector_index_.size();
    
    double base1 =  current_vector_*scale1;
    double base2 =  current_atomindex_*scale1*scale2;
    double base3 =  current_subvector_*scale1*scale2*scale3;
    return base1 + base2 + base3;
}
// end of file
