/** \file
This file contains a class which implements the scattering calculation for multipole moment based orientally averaged scattering ( which is all type scattering ).

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "scatter_devices/multipole_scatter_device.hpp"

// standard header
#include <algorithm>
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions.hpp>

// other headers
#include "math/coor3d.hpp"
#include "math/smath.hpp"
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "stager/data_stager.hpp"

using namespace std;

MPSphereScatterDevice::MPSphereScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<CartesianCoor3D> vectors,
    size_t NAF,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    AbstractScatterDevice(
        allcomm,
        partitioncomm,
        sample,
        vectors,
        NAF,
        fileservice_endpoint,
        monitorservice_endpoint
    ),
    current_moment_(0)
{
	sample_.coordinate_sets.set_representation(SPHERICAL);	
	
	fftw_complex* wspace= NULL;
    fftw_planF_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_planB_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void MPSphereScatterDevice::print_pre_stage_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Staging data...");
    }
}

void MPSphereScatterDevice::print_post_stage_info() { }

void MPSphereScatterDevice::print_pre_runner_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Starting computation...");
    }
}

void MPSphereScatterDevice::print_post_runner_info() { }

double MPSphereScatterDevice::progress() {
    double scale1 = 1.0/vectors_.size();
    double scale2 = 1.0/multipole_index_.size();
    
    double base1 =  current_vector_*scale1;
    double base2 =  current_moment_*scale1*scale2;
    return base1 + base2;
}

bool MPSphereScatterDevice::ram_check() {
	// inherit ram requirements for parent class
	bool state = AbstractScatterDevice::ram_check();
	
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
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.signal_buffer));
	Info::Inst()->write(string("required limits.computation.memory.exchange_buffer=")+
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.exchange_buffer));
	Info::Inst()->write(string("required limits.computation.memory.alignpad_buffer=")+
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.alignpad_buffer));
	}
	// final buffer for reduction:
	// NF*sizeof(fftw_complex)
	// however, at that moment the signal_buffer is deallocated
	// so we can ignore this requirement
	// since we're interested in the potential peak consumption..

	return state;
}

void MPSphereScatterDevice::init_moments(CartesianCoor3D& q) {
    multipole_index_.clear();	
    
    for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.multipole.moments.size(); ++i)
	{
		multipole_index_.push_back(Params::Inst()->scattering.average.orientation.multipole.moments[i]);
	}
    
    qvector_ = q;
    NM=multipole_index_.size();
}

void MPSphereScatterDevice::stage_data() {
    Timer& timer = timer_[boost::this_thread::get_id()];
    if (allcomm_.rank()==0) Info::Inst()->write(string("Forcing stager.mode=frames"));
    DataStagerByFrame data_stager(sample_,allcomm_,partitioncomm_,timer);
    p_coordinates = data_stager.stage();
}

MPSphereScatterDevice::~MPSphereScatterDevice() {
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

void MPSphereScatterDevice::worker() {
      
    workerbarrier->wait();

    Timer& timer = timer_[boost::this_thread::get_id()];
     
    while (true) {
        size_t moment_index;
        timer.start("sd:worker:wait");
        atscatter_.wait_and_pop(moment_index);  
        timer.stop("sd:worker:wait");
        
//        if (moment_index<NM) {
            timer.start("sd:worker:scatter");
            scatter(moment_index);
            timer.stop("sd:worker:scatter");
//        }

        workerbarrier->wait();
    }
}

fftw_complex* MPSphereScatterDevice::exchange() {
    
    size_t NNPP = partitioncomm_.size();
  
    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
  
    fftw_complex* atOUT = (fftw_complex*) fftw_malloc(NMAXF*NNPP*sizeof(fftw_complex)); 

    double* pIN = (double*) at_;
    double* pOUT = (double*) atOUT;

    boost::mpi::all_to_all(partitioncomm_,pIN,2*NMAXF,pOUT);

    return atOUT;
}

fftw_complex* MPSphereScatterDevice::alignpad(fftw_complex* at) {

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

void MPSphereScatterDevice::dsp(fftw_complex* at) {
    
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

void MPSphereScatterDevice::store(fftw_complex* at) {
    complex<double> a = smath::reduce<double>(at,NF) * (1.0/NF);
    afinal_ += a;
    a2final_ += a*conj(a);
    smath::add_elements(atfinal_,at,NF);
}

void MPSphereScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    Timer& timer = timer_[boost::this_thread::get_id()];

    timer.start("sd:c:init");
    init_moments(q);
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
    timer.stop("sd:c:init");

    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
    size_t NNPP = partitioncomm_.size();
            
    current_moment_=0;
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
            
            current_moment_+=NTHREADS;
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
            
            current_moment_+=std::min(NNPP,NM-i);
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
            p_atlocal = (double*) &(atlocal_[0][0]);;
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
    
    double factor = 1.0/(4*M_PI);
    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
        afinal_ *= factor;        
        a2final_ *= factor;                
    }
}

void MPSphereScatterDevice::scatterblock(size_t index,size_t count) {
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


void MPSphereScatterDevice::scatter(size_t this_moment) {
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
       offset = (this_moment%NTHREADS)*NF;
   } else {
       offset = (this_moment%NNPP)*NMAXF;       
   }

	DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);
   size_t NMYF = assignment.size();
   fftw_complex* p_a = &(at_[offset]);
   
   std::pair<long,long> moment = multipole_index_[this_moment];
   double ql = qvector_.length();
   double M_PI_four = 4*M_PI;

   SphericalCoor3D qs(qvector_);

   long l = moment.first;
   long m = moment.second;
   
   if (abs(m)>l) {
       Err::Inst()->write(string("Combination of Major and minor moment not allowed: l=") + boost::lexical_cast<string>(l) + string(", m")+boost::lexical_cast<string>(m));
       Err::Inst()->write("Consult the manual/reference.");
       throw;
   }
   
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       coor_t* p_data = &(p_coordinates[fi*NA*3]);

       complex<double> A = 0;
       for(size_t j = 0; j < NA; ++j) {
           double r = p_data[3*j];
           double phi = p_data[3*j+1];
           double theta = p_data[3*j+2];
       
           double esf = sfs[j];
           double p = ql*r;
//           cout << "r=" << r << endl;
//           cout << "phi=" << phi << endl;
//           cout << "theta=" << theta << endl;
   
	       complex<double> fmpiilesf = M_PI_four*pow(complex<double>(0,1.0),l) * esf;
	       complex<double> aabess = boost::math::sph_bessel(l,p);
//           cout << "fmpiilesf=" << fmpiilesf << endl;
//           cout << "aabess=" << aabess << endl;

	       complex<double> aa = conj(boost::math::spherical_harmonic(l,m,theta,phi)); 
//           cout << "aa=" << aa << endl;
           
           A += fmpiilesf * aabess* aa;// * boost::math::spherical_harmonic(l,m,qs.theta,qs.phi);
//           cout << "j=" << j << ", A=" << A << endl;
       }
       p_a[fi][0] = A.real();
       p_a[fi][1] = A.imag();
//       cout << "l,m=" << l << "," << m <<  ", A=" << A << endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
// Cylinder Multipole
////////////////////////////////////////////////////////////////////////////////

MPCylinderScatterDevice::MPCylinderScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<CartesianCoor3D> vectors,
    size_t NAF,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    AbstractScatterDevice(
        allcomm,
        partitioncomm,
        sample,
        vectors,
        NAF,
        fileservice_endpoint,
        monitorservice_endpoint
    ),
    current_moment_(0)
{
	sample_.coordinate_sets.set_representation(CYLINDRICAL);	
	
	fftw_complex* wspace= NULL;
    fftw_planF_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_planB_ = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void MPCylinderScatterDevice::print_pre_stage_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Staging data...");
    }
}

void MPCylinderScatterDevice::print_post_stage_info() { }

void MPCylinderScatterDevice::print_pre_runner_info() {
    if (allcomm_.rank()==0) {
        Info::Inst()->write("Starting computation...");
    }
}

void MPCylinderScatterDevice::print_post_runner_info() { }

double MPCylinderScatterDevice::progress() {
    double scale1 = 1.0/vectors_.size();
    double scale2 = 1.0/multipole_index_.size();
    
    double base1 =  current_vector_*scale1;
    double base2 =  current_moment_*scale1*scale2;
    return base1 + base2;
}


bool MPCylinderScatterDevice::ram_check() {
	// inherit ram requirements for parent class
	bool state = AbstractScatterDevice::ram_check();
	
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
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.signal_buffer));
	Info::Inst()->write(string("required limits.computation.memory.exchange_buffer=")+
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.exchange_buffer));
	Info::Inst()->write(string("required limits.computation.memory.alignpad_buffer=")+
		boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.alignpad_buffer));
	}
	// final buffer for reduction:
	// NF*sizeof(fftw_complex)
	// however, at that moment the signal_buffer is deallocated
	// so we can ignore this requirement
	// since we're interested in the potential peak consumption..

	return state;
}

void MPCylinderScatterDevice::init_moments(CartesianCoor3D& q) {
    multipole_index_.clear();	
    
    for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.multipole.moments.size(); ++i)
	{
		multipole_index_.push_back(Params::Inst()->scattering.average.orientation.multipole.moments[i]);
	}
    
    NM = multipole_index_.size();
    qvector_ = q;
}

void MPCylinderScatterDevice::stage_data() {
    Timer& timer = timer_[boost::this_thread::get_id()];
    DataStagerByFrame data_stager(sample_,allcomm_,partitioncomm_,timer);
    p_coordinates = data_stager.stage();
}

MPCylinderScatterDevice::~MPCylinderScatterDevice() {
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

void MPCylinderScatterDevice::worker() {
      
    workerbarrier->wait();

    Timer& timer = timer_[boost::this_thread::get_id()];
     
    while (true) {
        size_t moment_index;
        timer.start("sd:worker:wait");
        atscatter_.wait_and_pop(moment_index);  
        timer.stop("sd:worker:wait");
        
//        if (moment_index<NM) {
            timer.start("sd:worker:scatter");
            scatter(moment_index);
            timer.stop("sd:worker:scatter");
//        }

        workerbarrier->wait();
    }
}

fftw_complex* MPCylinderScatterDevice::exchange() {
    
    size_t NNPP = partitioncomm_.size();
  
    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
  
    fftw_complex* atOUT = (fftw_complex*) fftw_malloc(NMAXF*NNPP*sizeof(fftw_complex)); 

    double* pIN = (double*) at_;
    double* pOUT = (double*) atOUT;

    boost::mpi::all_to_all(partitioncomm_,pIN,2*NMAXF,pOUT);

    return atOUT;
}

fftw_complex* MPCylinderScatterDevice::alignpad(fftw_complex* at) {

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

void MPCylinderScatterDevice::dsp(fftw_complex* at) {
    
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

void MPCylinderScatterDevice::store(fftw_complex* at) {
    complex<double> a = smath::reduce<double>(at,NF) * (1.0/NF);
    afinal_ += a;
    a2final_ += a*conj(a);
    smath::add_elements(atfinal_,at,NF);
}

void MPCylinderScatterDevice::compute() {
    CartesianCoor3D q=vectors_[current_vector_];

    Timer& timer = timer_[boost::this_thread::get_id()];
    timer.start("sd:c:init");
    init_moments(q);
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
    timer.stop("sd:c:init");

    DivAssignment zeronode_assignment(partitioncomm_.size(),0,NF);
    size_t NMAXF = zeronode_assignment.max();
    size_t NNPP = partitioncomm_.size();
            
    current_moment_=0;
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
            
            current_moment_+=NTHREADS;
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
            
            current_moment_+=std::min(NNPP,NM-i);
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
    if (NN>1) {

        double* p_atfinal = (double*) &(atfinal_[0][0]);
        double* p_atlocal=NULL;
        if (partitioncomm_.rank()==0) {
            fftw_complex* atlocal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
            p_atlocal = (double*) &(atlocal_[0][0]);;
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
    
    double factor = 1.0/(2*M_PI);
    if (partitioncomm_.rank()==0) {
        smath::multiply_elements(factor,atfinal_,NF);
        afinal_ *= factor;        
        a2final_ *= factor;                
    }
}

void MPCylinderScatterDevice::scatterblock(size_t index,size_t count) {
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



void MPCylinderScatterDevice::scatter(size_t this_moment) {
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
       offset = (this_moment%NTHREADS)*NF;
   } else {
       offset = (this_moment%NNPP)*NMAXF;       
   }

   DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);
   size_t NMYF = assignment.size();
   fftw_complex* p_a = &(at_[offset]);
   
   std::pair<long,long> moment = multipole_index_[this_moment];

   CartesianCoor3D o = Params::Inst()->scattering.average.orientation.axis;
	
    // constructs a base out of thin air
   CartesianVectorBase base(o);
   CartesianCoor3D qprojected = base.project(qvector_);
   CylinderCoor3D qcylinder(qprojected);
   double qr = qcylinder.r;
   double qphi = qcylinder.phi;
   double qz = qcylinder.z;

   long l = moment.first;
   long m = moment.second;
   
   
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       coor_t* p_data = &(p_coordinates[fi*NA*3]);

       complex<double> A = 0;
       for(size_t j = 0; j < NA; ++j) {
           double r = p_data[3*j];
           double phi = p_data[3*j+1];
           double z = p_data[3*j+2];
       
           double esf = sfs[j];
           
           double parallel_sign = 1.0;
	       if ((z!=0) && (qz!=0)) {
	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
	       }

	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
	       double p = r*qr;
           double psiphi = phi-qphi; // review this!
           
           if (m==0 && l==0) {
               A += expi * (double)boost::math::cyl_bessel_j(0,p) * esf ;
           } else if (m==0) {    	       
    	       // if func==1
               complex<double> fac1 = 2.0*powf(-1.0,l)*boost::math::cyl_bessel_j(2*l,p);
               A += sqrt(0.5)*fac1 * expi *cos(2*l*psiphi) * esf;               
           } else if (m==1) {
               // if func==2
        	   complex<double> fac1 = 2.0*powf(-1.0,l)*boost::math::cyl_bessel_j(2*l,p);
               A += sqrt(0.5)*fac1 * expi *sin(2*l*psiphi) * esf;          
           } else if (m==2) {
               // if func==3
               complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*boost::math::cyl_bessel_j(2*l-1,p));
               A += sqrt(0.5)*fac2 * expi *cos((2*l-1)*psiphi) * esf;    
           } else if (m==3) {
               // if func==4
               complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*boost::math::cyl_bessel_j(2*l-1,p));
               A += sqrt(0.5)*fac2 * expi *sin((2*l-1)*psiphi) * esf;
           }
       }
       // normalize with integral over phi
       double norm = sqrt(2*M_PI);
       p_a[fi][0] = norm*A.real();
       p_a[fi][1] = norm*A.imag();
//       cout << "l,m=" << l << "," << m <<  ", A=" << A << endl;
   }
}

// end of file
