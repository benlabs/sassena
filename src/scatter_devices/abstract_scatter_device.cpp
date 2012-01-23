/** \file
This file contains the interface definition for all scattering devices and implements an abstract scattering device from which all other devices are derived.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
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
#include "exceptions/exceptions.hpp"
#include "math/coor3d.hpp"
#include "decomposition/assignment.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "stager/data_stager.hpp"

using namespace std;

AbstractScatterDevice::AbstractScatterDevice(
    boost::mpi::communicator allcomm,
    boost::mpi::communicator partitioncomm,
    Sample& sample,
    std::vector<CartesianCoor3D> vectors,
    size_t NAF,
    boost::asio::ip::tcp::endpoint fileservice_endpoint,
	boost::asio::ip::tcp::endpoint monitorservice_endpoint
) :
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),
    sample_(sample),
    vectors_(vectors),
    current_vector_(0),
    atfinal_(NULL),
    afinal_(0),
    a2final_(0)
{
    
    p_hdf5writer_ = boost::shared_ptr<HDF5WriterClient>(new HDF5WriterClient(fileservice_endpoint));
    p_monitor_ = boost::shared_ptr<MonitorClient>(new MonitorClient(monitorservice_endpoint));
    	
    NN = partitioncomm_.size();
	NF = sample_.coordinate_sets.size();
	std::string target = Params::Inst()->stager.target;
	NA = sample_.atoms.selections[target]->size(); // Number of Atoms
	
	sample_.coordinate_sets.set_selection(sample.atoms.selections[target]);
	
	scatterfactors.set_sample(sample_);
	scatterfactors.set_selection(sample_.atoms.selections[target]);
	scatterfactors.set_background(true);
	
	// defer memory allocation for scattering data
    atfinal_ = NULL;
    
    Timer blank_timer;
    timer_.insert(map<boost::thread::id,Timer>::value_type(boost::this_thread::get_id(),blank_timer));	    
}

AbstractScatterDevice::~AbstractScatterDevice() {
    if (atfinal_!=NULL) {
        fftw_free(atfinal_);
        atfinal_=NULL;
    }
}

bool AbstractScatterDevice::ram_check() {
	bool state = true;
	
	size_t memscale = Params::Inst()->limits.computation.memory.scale;
	
	// atfinal_
	size_t bytesize_result_buffer=NF*sizeof(fftw_complex);
	
    if (bytesize_result_buffer>memscale*Params::Inst()->limits.computation.memory.result_buffer) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("limits.computation.memory.result_buffer too small");
			Err::Inst()->write(string("limits.computation.memory.result_buffer=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.result_buffer));
			Err::Inst()->write(string("limits.computation.memory.scale=")+ boost::lexical_cast<string>(Params::Inst()->limits.computation.memory.scale));	
            Err::Inst()->write(string("requested: ")+boost::lexical_cast<string>(bytesize_result_buffer));            
        }
		state = false;
    }
	if (allcomm_.rank()==0) {
	Info::Inst()->write(string("required limits.computation.memory.result_buffer=")+
		boost::lexical_cast<string>(bytesize_result_buffer));
	}
	return state;
}

void AbstractScatterDevice::run() {
    
    Timer& timer = timer_[boost::this_thread::get_id()];

	// check memory requirements here
	if (allcomm_.rank()==0) Info::Inst()->write("Checking memory requirements for scattering calculation...");
	if (!ram_check()) throw sassena::terminate_request();	
	if (allcomm_.rank()==0) Info::Inst()->write("Memory limits ok");

    print_pre_stage_info();
    allcomm_.barrier();
    timer.start("sd:stage");
    stage_data();
    timer.stop("sd:stage");
    allcomm_.barrier();    
    print_post_stage_info();

    if (allcomm_.rank()==0) {
        scatterfactors.update(CartesianCoor3D(0,0,0));
        Info::Inst()->write("Target initialized. ");
        Info::Inst()->write(string("Target produces a background scattering length density of ")+boost::lexical_cast<string>(scatterfactors.compute_background(CartesianCoor3D(0,0,0))));
    }

    // allocate memory for computation now.
    atfinal_ = (fftw_complex*) fftw_malloc(NF*sizeof(fftw_complex));
    memset(atfinal_,0,NF*sizeof(fftw_complex));

    print_pre_runner_info();
    allcomm_.barrier();
    timer.start("sd:runner");
    runner();
    timer.stop("sd:runner");
    allcomm_.barrier();
    print_post_runner_info();

    // notify the services to finalize
    p_monitor_->update(allcomm_.rank(),1.0);
    p_hdf5writer_->flush();        
}

void AbstractScatterDevice::runner() {

    Timer& timer = timer_[boost::this_thread::get_id()];
 
     start_workers();
     
     size_t samplingfactor = Params::Inst()->limits.services.monitor.sampling;
     if (samplingfactor==0) {
         samplingfactor = allcomm_.size()/20;
         if (allcomm_.size()<20) samplingfactor=1;
     }
     p_monitor_->set_samplingfactor(samplingfactor);
     if (allcomm_.rank()==0) {
         Info::Inst()->write(string("Setting progress sampling factor for monitoring to ")+boost::lexical_cast<string>(samplingfactor));
         p_monitor_->set_samplingfactor_server(samplingfactor);
         p_monitor_->reset_server();
     }
     while(status()==0) {
            
         timer.start("sd:compute");
 	     compute();
         timer.stop("sd:compute"); 
            
         timer.start("sd:write");        	
	     write();
	     timer.stop("sd:write");
            		    
        next();
     }             
     stop_workers();
}


void AbstractScatterDevice::start_workers() {

    size_t numthreads = Params::Inst()->limits.computation.threads;
//    if (worker_threads>NM) {
//        if (allcomm_.rank()==0) {
//            Warn::Inst()->write("Number of threads limited by resolution for averaging");
//            Warn::Inst()->write(string("Number of moments (resolution): ")+boost::lexical_cast<string>(NM));
//            Warn::Inst()->write(string("Number of threads requested: ")+boost::lexical_cast<string>(worker_threads));
//            Warn::Inst()->write(string("Setting limits.computation.threads.worker=")+boost::lexical_cast<string>(NM));
//        }
//        worker_threads = NM;
//    }
    
    if ( (partitioncomm_.size()!=1) &&  (partitioncomm_.size()<numthreads) ) {
        if (allcomm_.rank()==0) {
            Warn::Inst()->write("Partition smaller than number of threads.");
            Warn::Inst()->write("Can not utilize more threads than partition size.");
            Warn::Inst()->write(string("Setting limits.computation.threads.worker=")+boost::lexical_cast<string>(partitioncomm_.size()));
        }
        numthreads = partitioncomm_.size();
    }

    workerbarrier = new boost::barrier(numthreads+1);

    Timer blank_timer;
    for(size_t i = 0; i < numthreads; ++i)
    {
        boost::thread* p_t = new boost::thread(boost::bind(&AbstractScatterDevice::worker,this));
        worker_threads.push(p_t);
        timer_.insert(map<boost::thread::id,Timer>::value_type(p_t->get_id(),blank_timer));	
    }
    
    workerbarrier->wait();
}

void AbstractScatterDevice::stop_workers() {

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
    delete workerbarrier;
}

void AbstractScatterDevice::next() {
	if (current_vector_>=vectors_.size()) return;
	current_vector_++;    
}

double AbstractScatterDevice::progress() {
    double scale = 1.0/vectors_.size();
    double base =  current_vector_*scale;
    return base;    
}

size_t AbstractScatterDevice::status() {
    if (current_vector_==vectors_.size()) return 1; else return 0;
}

void AbstractScatterDevice::write() {
    if (partitioncomm_.rank()==0) {
        CartesianCoor3D vector = vectors_[current_vector_]; 
        p_hdf5writer_->write(vector,atfinal_,NF,afinal_,a2final_);
    }
}

std::map<boost::thread::id,Timer>& AbstractScatterDevice::getTimer() {
    // collect timer information from allcomm_
    
    return timer_;    
}

// end of file
