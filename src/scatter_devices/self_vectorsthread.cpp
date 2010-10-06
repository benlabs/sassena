/*
 *  self_vectorsthread.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/self_vectorsthread.hpp"

// standard header
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
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "scatter_devices/data_stager.hpp"

using namespace std;

// non-specific within VECTORS
SelfVectorsThreadScatterDevice::SelfVectorsThreadScatterDevice(
		boost::mpi::communicator scatter_comm,
		boost::mpi::communicator fqt_comm,
		Sample& sample,
		vector<pair<size_t,CartesianCoor3D> > QIV,
		vector<size_t> assignment,
		boost::asio::ip::tcp::endpoint fileservice_endpoint,
		boost::asio::ip::tcp::endpoint monitorservice_endpoint)
	: p_coordinates(NULL)
{
	m_qvectorindexpairs= QIV;
	m_current_qvector =0;

	m_scattercomm = scatter_comm;
	m_fqtcomm = fqt_comm;
    
    m_assignment = assignment;
	p_sample = &sample; // keep a reference to the sample

    p_hdf5writer = new HDF5WriterClient(fileservice_endpoint);
    p_monitor = new MonitorClient(monitorservice_endpoint);

	string target = Params::Inst()->scattering.target;

	NN = m_fqtcomm.size(); // Number of Nodes
	NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	NF = sample.coordinate_sets.size();
    NM = Params::Inst()->scattering.average.orientation.vectors.size();
		
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

// SPECIFIC
void SelfVectorsThreadScatterDevice::stage_data() {
    DataStagerByAtom data_stager(*p_sample,m_scattercomm,m_assignment);
    p_coordinates = data_stager.stage();
}

// non-specific
SelfVectorsThreadScatterDevice::~SelfVectorsThreadScatterDevice() {
    if (p_coordinates!=NULL) free(p_coordinates);
}

// non-specific
double SelfVectorsThreadScatterDevice::progress() {
    double scale1 = 1.0/m_qvectorindexpairs.size();
    double scale2 = 1.0/(NM*m_assignment.size());
    
    double base =  m_current_qvector*scale1;
    double rel = worker2_counter*scale2;
    return base + scale1*rel;    
}

// non-specific
size_t SelfVectorsThreadScatterDevice::status() {
    return m_current_qvector/m_qvectorindexpairs.size();
}

// non-specific
void SelfVectorsThreadScatterDevice::run() {
    
    if (m_scattercomm.rank()==0) {
        Info::Inst()->write("Staging data...");
    }
    
    stage_data();
    	
    
    if (Params::Inst()->debug.print.orientations) {     
    if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        if (m_scattercomm.rank()==0) {		
    		Info::Inst()->write("Qvectors orientations used for averaging: ");
            for (size_t i=0;i<Params::Inst()->scattering.average.orientation.vectors.size();i++) {
                string qvector = "";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].x);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].y);
                qvector += " ";
                qvector += boost::lexical_cast<string>(Params::Inst()->scattering.average.orientation.vectors[i].z);
    		    Info::Inst()->write(qvector);                                
            }
    	}	    
    }    
    }
    
    if (m_scattercomm.rank()==0) {
        Info::Inst()->write("Starting computation...");
    }
    
    runner();
    // notify the services to finalize
    if (Params::Inst()->debug.monitor.progress) {
        p_monitor->update(m_scattercomm.rank(),1.0);
    }
    if (Params::Inst()->debug.iowrite) {
        p_hdf5writer->flush();        
    }
}

// non-specific
void SelfVectorsThreadScatterDevice::runner() {
 
    bool threads_on = Params::Inst()->limits.computation.threads;
    
    if (threads_on) start_workers();

    while(status()==0) {
    	compute(!threads_on);
    	if (Params::Inst()->debug.iowrite) {
    		write();		    
    	}
        next();
    }    

    if (threads_on) stop_workers(); 
}

// non-specific, IF worker1,worker2,worker3 is there
void SelfVectorsThreadScatterDevice::start_workers() {

    size_t worker2_threads = 1;
    size_t worker1_threads = Params::Inst()->limits.computation.worker1_threads;
    if (worker1_threads==0) {
        size_t busy = worker2_threads;
        if (busy>=boost::thread::hardware_concurrency()) {
            worker1_threads = 1;
        } else {
            worker1_threads = boost::thread::hardware_concurrency() - busy;
        }
    }
    
    // don't allow more worker1 threads than NM are available!
    if (worker1_threads>NM) worker1_threads = NM;
    
    
    for(size_t i = 0; i < worker1_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&SelfVectorsThreadScatterDevice::worker1,this,true) ));	
    }

    for(size_t i = 0; i < worker2_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&SelfVectorsThreadScatterDevice::worker2,this,true) ));	
    }
}

// non-specific
void SelfVectorsThreadScatterDevice::stop_workers() {

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
}


// non-specific
void SelfVectorsThreadScatterDevice::next() {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	m_current_qvector+=1;    
}

// non-specific
void SelfVectorsThreadScatterDevice::write() {
    if (Params::Inst()->debug.iowrite) {
        if ((m_fqtcomm.rank()==0)) {
            size_t qindex = m_qvectorindexpairs[m_current_qvector].first; 
            p_hdf5writer->write(qindex,m_spectrum);
        }        
    }
}

// SPECIFIC
void SelfVectorsThreadScatterDevice::compute(bool marshal) {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	CartesianCoor3D q=m_qvectorindexpairs[m_current_qvector].second;

    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
            
    m_spectrum.clear();
    atfinal.clear();
    m_spectrum.resize(NF,0);
    atfinal.resize(NF,0);
    
    worker2_counter = 0;
    worker2_done = false;
    boost::mutex::scoped_lock w2l(worker2_mutex);

    if (marshal) {
        for(size_t i = 0; i < NM; ++i)
        {
            for(size_t n = 0; n < m_assignment.size(); ++n)
            {
                at0.push(make_pair(i,n));
                worker1(false);
                worker2(false);                                   
            }
        }
        if (Params::Inst()->debug.monitor.progress) {
            p_monitor->update(m_scattercomm.rank(),progress());    						            
        }
    } else {
        for(size_t i = 0; i < NM; ++i) {
            for(size_t n = 0; n < m_assignment.size(); ++n)
            {
                at0.push(make_pair(i,n));
            }
        }
        while (!worker2_done)  worker2_notifier.wait(w2l);
    }   

    
    vector<complex<double> > at_local(NF,0);
    
    double* p_atfinal = (double*) &(atfinal[0]);
    double* p_atlocal = (double*) &(at_local[0]);
    
    if (NN>1) {
        boost::mpi::reduce(m_fqtcomm,p_atfinal,2*NF,p_atlocal,std::plus<double>(),0);
        for(size_t j = 0; j < NF; ++j)
        {
            m_spectrum[j] += std::complex<double>((at_local[j].real())/NM,(at_local[j].imag())/NM);
        }         			
    } 
    else {
        for(size_t j = 0; j < NF; ++j)
        {
            m_spectrum[j] += std::complex<double>((atfinal[j].real())/NM,(atfinal[j].imag())/NM);
        }         			
    }
    

}

// non-specific
void SelfVectorsThreadScatterDevice::init(CartesianCoor3D& q) {
    qvectors.clear();	
	
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
		if (Params::Inst()->scattering.average.orientation.vectors.type=="sphere") {
			double ql = q.length();
			
			for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
			{
				qvectors.push_back(ql*Params::Inst()->scattering.average.orientation.vectors[i]);
			}					
		} else if (Params::Inst()->scattering.average.orientation.vectors.type=="cylinder") {
			CartesianCoor3D o = Params::Inst()->scattering.average.orientation.vectors.axis;
			// make sure o is normalized;
			o = o / o.length();
			
			// get the part of the scattering vector perpenticular to the o- orientation
			CartesianCoor3D qparallel = (o*q)*o; 
			CartesianCoor3D qperpenticular = q - qparallel; 			
			double qperpenticular_l = qperpenticular.length();
			//double qparallel_l = qparallel.length();

			CartesianCoor3D e1 = o.cross_product(qperpenticular) ;
			CartesianCoor3D e2 = qperpenticular;
			
			if (qperpenticular_l==0.0)  { 
				qvectors.push_back( qparallel );
			}
			else {
				for(size_t i = 0; i < Params::Inst()->scattering.average.orientation.vectors.size(); ++i)
				{
					CartesianCoor3D& vec = Params::Inst()->scattering.average.orientation.vectors[i];
					CartesianCoor3D qnew = qparallel + vec.x * e1 + vec.y * e2;					
					qvectors.push_back(qnew);
				}
			}				
		}
	} else {
		qvectors.push_back(q);
	}
}

////////////////////////////////////////////////////////////////////////////////
//  specific functions
////////////////////////////////////////////////////////////////////////////////

void SelfVectorsThreadScatterDevice::worker1(bool loop) {
    
    size_t maxbuffer_bytesize = Params::Inst()->limits.memory.at1_buffer;
    size_t timeout = Params::Inst()->limits.times.at1_buffer;

    size_t single_entry_bytesize = NF*2*sizeof(double);
    size_t max_entries = maxbuffer_bytesize/single_entry_bytesize;
    
    while (true) {
        
        // retrieve from queue
        std::pair<size_t,size_t> qmoment_atom;

        // QoS here
        while (at1.size()>max_entries) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(timeout));
        }
        at0.wait_and_pop(qmoment_atom);
        // operation
        std::vector<std::complex<double> >* p_at_local = scatter(qmoment_atom.first,qmoment_atom.second);
        
        // apply dsp
	    if (Params::Inst()->scattering.dsp.type=="autocorrelate") {
            if (Params::Inst()->scattering.dsp.method=="direct") {
                smath::auto_correlate_direct(*p_at_local);	        
            } else if (Params::Inst()->scattering.dsp.method=="fftw") {
                smath::auto_correlate_fftw(*p_at_local);	        
            } else {
                Err::Inst()->write("Correlation method not understood");
                Err::Inst()->write("scattering.dsp.method == direct, fftw");        
                throw;
            }
		} else if (Params::Inst()->scattering.dsp.type=="square") {
            smath::square_elements(*p_at_local);
		}
		
        at1.push(p_at_local);
        if (!loop) {
            at1.push(NULL);
            break;
        }
    }
}

void SelfVectorsThreadScatterDevice::worker2(bool loop) {
    while ( true ) {
        std::vector< std::complex<double> >* p_at_local=NULL;
        if (loop) {
            at1.wait_and_pop(p_at_local);
            if (p_at_local==NULL) { // ats as sentinal
                worker2_done = true;
                boost::mutex::scoped_lock w2l(worker2_mutex);
                w2l.unlock();
                worker2_notifier.notify_all();
                       
                if (Params::Inst()->debug.monitor.progress) {
                   p_monitor->update(m_scattercomm.rank(),progress());    						            
                }         
                continue;
            }
        } else {
            at1.try_pop(p_at_local);            
            if(p_at_local==NULL) break;
        }
        worker2_counter++;
        smath::add_elements(atfinal,*p_at_local);
        if (Params::Inst()->debug.monitor.progress) {
           p_monitor->update(m_scattercomm.rank(),progress());    						            
        }         

        delete p_at_local;
    }
}

// computation


std::vector<std::complex<double> >* SelfVectorsThreadScatterDevice::scatter(size_t mi,size_t ai) {

	// this is broken <-- revise this!!!
	double s = scatterfactors.get(m_assignment[ai]);
	
    std::vector<std::complex<double> >* p_at_local = new std::vector<std::complex<double> >(NF);

    double qx = qvectors[mi].x;
    double qy = qvectors[mi].y;
    double qz = qvectors[mi].z;

    
    coor_t* p_data = &(p_coordinates[ai*NF*3]);

	for(size_t j = 0; j < NF; ++j)
	{
		coor_t x1 = p_data[j*3];
		coor_t y1 = p_data[j*3 + 1];
		coor_t z1 = p_data[j*3 + 2];

		double p1 = x1*qx+y1*qy+z1*qz;
        
        double sp1 = sin(p1);
        double cp1 = cos(p1);
        (*p_at_local)[j] = s*complex<double>(cp1,sp1);
	}
  
    return p_at_local;
}


// end of file
