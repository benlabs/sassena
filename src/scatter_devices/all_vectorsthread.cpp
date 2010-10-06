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
#include "scatter_devices/all_vectorsthread.hpp"

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
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"
#include "scatter_devices/data_stager.hpp"

using namespace std;

AllVectorsThreadScatterDevice::AllVectorsThreadScatterDevice(
		boost::mpi::communicator scatter_comm,
		boost::mpi::communicator fqt_comm,
		Sample& sample,
		std::vector<std::pair<size_t,CartesianCoor3D> > QIV,
        std::vector<size_t> assignment,
		boost::asio::ip::tcp::endpoint fileservice_endpoint,
		boost::asio::ip::tcp::endpoint monitorservice_endpoint)
	: p_coordinates(NULL)
{
	m_qvectorindexpairs= QIV;
	m_current_qvector =0;

	m_scattercomm = scatter_comm;
	m_fqtcomm = fqt_comm;
    
	p_sample = &sample; // keep a reference to the sample

    p_hdf5writer = new HDF5WriterClient(fileservice_endpoint);
    p_monitor = new MonitorClient(monitorservice_endpoint);

	string target = Params::Inst()->scattering.target;

	NN = m_fqtcomm.size(); // Number of Nodes
	NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	NF = sample.coordinate_sets.size();
    
    NM = Params::Inst()->scattering.average.orientation.vectors.size();
	
	size_t rank = m_fqtcomm.rank();

	EvenDecompose edecomp(NF,NN);
	
	myframes = edecomp.indexes_for(rank);
    // blocking factor:
    // max(nn,nq)
	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

    // overwrite representation style in subclass !!
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	//if (Params::Inst()->scattering.center) {
    //	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	//}

	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

void AllVectorsThreadScatterDevice::stage_data() {
    DataStagerByFrame data_stager(*p_sample,m_scattercomm,myframes);
    p_coordinates = data_stager.stage();
}

AllVectorsThreadScatterDevice::~AllVectorsThreadScatterDevice() {
    if (p_coordinates!=NULL) free(p_coordinates);
}

// scatter_frame implemented by concrete class
// scatter_frames implemented by concrete class

// init_averaging implemented by concrete class
// get_numberofmoments implemented by concrete class

// interface functions

void AllVectorsThreadScatterDevice::run() {
    
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

void AllVectorsThreadScatterDevice::runner() {
 
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


void AllVectorsThreadScatterDevice::start_workers() {

    size_t worker3_threads = 1; // fix this at the moment;
    size_t worker2_threads = 1;
    size_t worker1_threads = Params::Inst()->limits.computation.worker1_threads;
    if (worker1_threads==0) {
        size_t busy = worker2_threads+worker3_threads;
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
        worker_threads.push( new boost::thread(boost::bind(&AllVectorsThreadScatterDevice::worker1,this,true) ));	
    }

    worker_threads.push( new boost::thread(boost::bind(&AllVectorsThreadScatterDevice::worker2,this,true) ));	

    for(size_t i = 0; i < worker3_threads; ++i)
    {
        worker_threads.push( new boost::thread(boost::bind(&AllVectorsThreadScatterDevice::worker3,this,true) ));	
    }
}

void AllVectorsThreadScatterDevice::stop_workers() {

    while (worker_threads.size()>0) {
        boost::thread* p_thread = worker_threads.front();
        worker_threads.pop();
        p_thread->interrupt();
        delete p_thread;
    }
}


void AllVectorsThreadScatterDevice::worker1(bool loop) {
    
    size_t maxbuffer_bytesize = Params::Inst()->limits.memory.at1_buffer;
    size_t timeout = Params::Inst()->limits.times.at1_buffer;

    size_t single_entry_bytesize = myframes.size()*2*sizeof(double);
    size_t max_entries = maxbuffer_bytesize/single_entry_bytesize;
    
    while (true) {
        
        // retrieve from queue
        size_t qmoment;

        // QoS here
        while (at1.size()>max_entries) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(timeout));
        }
        at0.wait_and_pop(qmoment);  
        // operation
        scatter(qmoment);
        
        if (!loop) break;
    }
}


void AllVectorsThreadScatterDevice::worker2(bool loop) {
    
    size_t NN = m_fqtcomm.size();
    size_t NF = p_sample->coordinate_sets.size();
    size_t NMYF = myframes.size();    
    size_t MYRANK = m_fqtcomm.rank();
    
    // local data holders:
    EvenDecompose frameindexes(NF,NN);
    vector<size_t> nframe_offsets(NN);
    vector<size_t> nframes(NN);
    size_t nframe_integral=0;
    for (size_t i=0;i<NN;i++) {
        nframes[i]=frameindexes.indexes_for(i).size();
        nframe_offsets[i]=nframe_integral;
        nframe_integral+=nframes[i];        
    }
    
    std::map<size_t,std::vector<std::complex<double> >* > local_queue;
    
    while (true) {
        // retrieve from queue
        std::pair<size_t,std::vector< std::complex<double> >* > at_local_pair(0,NULL);
        std::vector< std::complex<double>  >* p_at_local = NULL;

        // first check local queue
        size_t key = worker2_counter;
        if (local_queue.find(key)==local_queue.end()) {
            at1.wait_and_pop(at_local_pair);   
        } else {
            at_local_pair.first = key;
            at_local_pair.second = local_queue[key];
            local_queue.erase(key);
        }
        
        size_t qindex = at_local_pair.first;
        p_at_local = at_local_pair.second;

        if (qindex!=key) {
            local_queue[at_local_pair.first]=at_local_pair.second;
            continue;
        }
        // operation
        
        // target node: vector_count % NN
        size_t target_node = (qindex%NN);
            
        if (MYRANK==target_node) {
            std::vector< std::complex<double> >* p_at_allframes = new std::vector<std::complex<double> >(NF);
                    
            for(size_t i = 0; i < NN; ++i)
            {
                double* pp_at_allframes = (double*) &((*p_at_allframes)[nframe_offsets[i]]);
                if (i==MYRANK) {
                    memcpy(pp_at_allframes,&((*p_at_local)[0]),2*nframes[i]*sizeof(double));
                } else {
                    m_fqtcomm.recv(i,0,pp_at_allframes,2*nframes[i]); 
                }
            }        
            
            at2.push(p_at_allframes);  
        } else {
            double* p_at_frames = (double*) &((*p_at_local)[0]);
            m_fqtcomm.send(target_node,0,p_at_frames,2*NMYF);
        }
        delete p_at_local;

        worker2_counter++;
        if (loop && (worker2_counter==NM)) {
            at2.push(NULL);
        }            
        
        if (!loop) break;
    }
    
}

void AllVectorsThreadScatterDevice::worker3(bool loop) {
    
    size_t counter=0;
    while ( true ) {
        std::vector< std::complex<double> >* p_at_local=NULL;
        if (loop) {
            at2.wait_and_pop(p_at_local);
            if (p_at_local==NULL) { // ats as sentinal
                worker3_done = true;
                counter=0;
                boost::mutex::scoped_lock w3l(worker3_mutex);
                w3l.unlock();
                worker3_notifier.notify_all();
                       
                if (Params::Inst()->debug.monitor.progress) {
                   double p = progress() + (1.0/m_qvectorindexpairs.size());
                   p_monitor->update(m_scattercomm.rank(),p);    						            
                }         
                continue;
            }
        } else {
            at2.try_pop(p_at_local);
            worker3_done = true;
            if (p_at_local==NULL) break;
        }
        
		// correlate or sum up
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
		} else if (Params::Inst()->scattering.dsp.type=="square")  {
            smath::square_elements(*p_at_local);
		}
		
        boost::mutex::scoped_lock at3l(at3_mutex);
        for(size_t n=0;n<at3.size();n++) {
            at3[n]+=(*p_at_local)[n];
        }
        at3l.unlock();

        delete p_at_local;
        counter++;
        
        if (Params::Inst()->debug.monitor.progress) {
            double p = progress() + (1.0/m_qvectorindexpairs.size())*(counter*1.0/NM);
            p_monitor->update(m_scattercomm.rank(),p);    						            
        }
                
        if (!loop) break;
	}	    
	
}

void AllVectorsThreadScatterDevice::next() {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	m_current_qvector+=1;    
}

double AllVectorsThreadScatterDevice::progress() {
    return m_current_qvector*1.0/m_qvectorindexpairs.size();    
}

size_t AllVectorsThreadScatterDevice::status() {
    return m_current_qvector/m_qvectorindexpairs.size();
}

void AllVectorsThreadScatterDevice::write() {
    if (Params::Inst()->debug.iowrite) {
        if ((m_fqtcomm.rank()==0)) {
            size_t qindex = m_qvectorindexpairs[m_current_qvector].first; 
            p_hdf5writer->write(qindex,m_spectrum);
        }        
    }
}

void AllVectorsThreadScatterDevice::compute(bool marshal) {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	CartesianCoor3D q=m_qvectorindexpairs[m_current_qvector].second;

    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
            
    m_spectrum.clear();
    at3.clear();
    m_spectrum.resize(NF,0);
    at3.resize(NF,0);
    
    worker2_counter = 0;
    worker3_done = false;
    boost::mutex::scoped_lock w3l(worker3_mutex);

    if (marshal) {    	
        for(size_t i = 0; i < NM; ++i)
        {
            at0.push(i);
            worker1(false);
            worker2(false);
            worker3(false);                       
        }
        if (Params::Inst()->debug.monitor.progress) {
            double p = progress();
            p_monitor->update(m_scattercomm.rank(),p);    						            
        }
        
    } else {
        for(size_t i = 0; i < NM; ++i) at0.push(i);
        while (!worker3_done)  worker3_notifier.wait(w3l);
    }

        

    
    vector<complex<double> > at_local(NF,0);
    
    double* p_at3 = (double*) &(at3[0]);
    double* p_atlocal = (double*) &(at_local[0]);
    
    if (NN>1) {
        boost::mpi::reduce(m_fqtcomm,p_at3,2*NF,p_atlocal,std::plus<double>(),0);
        for(size_t j = 0; j < NF; ++j)
        {
            m_spectrum[j] += std::complex<double>((at_local[j].real())/NM,(at_local[j].imag())/NM);
        }         			
    } 
    else {
        for(size_t j = 0; j < NF; ++j)
        {
            m_spectrum[j] += std::complex<double>((at3[j].real())/NM,(at3[j].imag())/NM);
        }         			
    }

}


void AllVectorsThreadScatterDevice::scatter(size_t moffset) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;
   
   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = myframes.size();

   std::vector<complex<double> >* p_a = new std::vector<complex<double> >(NMYF,0);
   
   string target = Params::Inst()->scattering.target;
   size_t NA = p_sample->atoms.selections[target].indexes.size();

   CartesianCoor3D& q = qvectors[moffset];
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       coor_t* p_data = &(p_coordinates[fi*NA*3]);
               
       double Ar =0; 
       double Ai = 0;
       coor_t x,y,z;
       double qx = q.x;
       double qy = q.y;
       double qz = q.z;
       double esf,p;
    
	   for(size_t j = 0; j < NA; ++j) {   		

           
	       esf = sfs[j];
	       x = p_data[3*j];
 	       y = p_data[3*j+1];
	       z = p_data[3*j+2];
       
           p =  x*qx + y*qy + z*qz;
           
           Ar += sfs[j]*cos(p);
           Ai += sfs[j]*sin(p);
	   }
       (*p_a)[fi] = complex<double>(Ar,Ai);
       //(*p_a)[fi] = complex<double>(2,3);
	}
    
    at1.push(make_pair(moffset,p_a));
}

void AllVectorsThreadScatterDevice::init(CartesianCoor3D& q) {
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
