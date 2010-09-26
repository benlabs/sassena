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
#include "threadpool/threadpool.hpp"

// other headers
#include "math/coor3d.hpp"
#include "math/smath.hpp"
#include "decomposition/decompose.hpp"
#include "scatter_devices/io/h5_fqt_interface.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;

AllVectorsThreadScatterDevice::AllVectorsThreadScatterDevice(
		boost::mpi::communicator scatter_comm,
		boost::mpi::communicator fqt_comm,
		Sample& sample,
		vector<pair<size_t,CartesianCoor3D> > QIV,
		string fqt_filename)
{
	m_qvectorindexpairs= QIV;
	m_current_qvector =0;
	m_fqt_filename = fqt_filename;

	m_scattercomm = scatter_comm;
	m_fqtcomm = fqt_comm;
    m_writeflag = false;

	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t NN = m_fqtcomm.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	size_t rank = m_fqtcomm.rank();

	EvenDecompose edecomp(NF,NN);
	
	myframes = edecomp.indexes_for(rank);
    // blocking factor:
    // max(nn,nq)

	if (m_fqtcomm.rank()==0) {

        size_t memusage_scatmat = 2*sizeof(double)*myframes.size()*1;
                
        size_t memusage_per_cs = 3*sizeof(double)*NA;
        size_t memusage_allcs = memusage_per_cs*myframes.size();
        
        
        Info::Inst()->write(string("Memory Requirements per node: "));
		Info::Inst()->write(string("Scattering Matrix: ")+to_s(memusage_scatmat)+string(" bytes"));


        // fault if not enough memory for scattering matrix
        if (memusage_scatmat>Params::Inst()->limits.memory.scattering_matrix) {
			Err::Inst()->write(string("Insufficient Buffer size for scattering matrix."));            
			Err::Inst()->write(string("Size required:")+to_s(memusage_scatmat)+string(" bytes"));            
			Err::Inst()->write(string("Configuration Parameter: limits.memory.scattering_matrix"));
            throw;
        }

		Info::Inst()->write(string("Coordinate Sets: ")+to_s(myframes.size()*memusage_per_cs)+string(" bytes"));		

        // warn if not enough memory for coordinate sets (cacheable system)
		if (Params::Inst()->runtime.limits.cache.coordinate_sets<myframes.size()) {
			Err::Inst()->write(string("Insufficient Buffer size for coordinate sets."));
			Err::Inst()->write(string("Need at least: ")+to_s(memusage_allcs)+string(" bytes"));
			Err::Inst()->write(string("Configuration Parameter: limits.memory.coordinate_sets"));
            throw;
		}		
	}
	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

    // overwrite representation style in subclass !!
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	//if (Params::Inst()->scattering.center) {
    //	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	//}
	
    for (size_t mf=0;mf<myframes.size();mf++) {
        csets.push_back(p_sample->coordinate_sets.load(myframes[mf]));
    }

	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
	
    size_t NMYF = myframes.size();
	
    size_t NM = 1;
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        NM = Params::Inst()->scattering.average.orientation.vectors.size();
    }
    size_t NMBLOCK = (NN<NM) ? NN : NM;
    
	if (m_fqtcomm.rank()==0) {

        size_t memusage_scatmat = 2*sizeof(double)*NMBLOCK*NMYF;

		Info::Inst()->write(string("Memory(Scattering Matrix): ")+to_s(memusage_scatmat)+string(" bytes"));

        // fault if not enough memory for scattering matrix
        if (memusage_scatmat>Params::Inst()->limits.memory.scattering_matrix) {
			Err::Inst()->write(string("Insufficient Buffer size for scattering matrix."));            
			Err::Inst()->write(string("Size required:")+to_s(memusage_scatmat)+string(" bytes"));            
			Err::Inst()->write(string("Configuration Parameter: limits.memory.scattering_matrix"));            
        }
	}
}

AllVectorsThreadScatterDevice::~AllVectorsThreadScatterDevice() {

}

// scatter_frame implemented by concrete class
// scatter_frames implemented by concrete class

// init_averaging implemented by concrete class
// get_numberofmoments implemented by concrete class

void AllVectorsThreadScatterDevice::correlate(std::vector<std::complex<double> >& data) {

    size_t NF = data.size();
    
    // make a local copy to allow to override data
    std::vector<std::complex<double> > data_local = data;
    
    std::vector< std::complex<double> >& complete_a = data_local;
    std::vector< std::complex<double> >& correlated_a = data;
    
    if (Params::Inst()->scattering.correlation.method=="direct") {
        
        // direct
        for(size_t tau = 0; tau < NF; ++tau)
        {
            correlated_a[tau] = 0;
        	size_t last_starting_frame = NF-tau;
        	for(size_t k = 0; k < last_starting_frame; ++k)
        	{
        		complex<double>& a1 = complete_a[k];
        		complex<double>& a2 = complete_a[k+tau];
        		correlated_a[tau] += conj(a1)*a2;
        	}
        	correlated_a[tau] /= (last_starting_frame); 		
        }
        
    } else if (Params::Inst()->scattering.correlation.method=="fftw") {
        
        fftw_complex *wspace;
        fftw_plan p1,p2;
        wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*NF);
        p1 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
        p2 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_BACKWARD, FFTW_ESTIMATE);


        for(size_t i = 0; i < NF; ++i) {
            wspace[i][0]= complete_a[i].real();            
            wspace[i][1]= complete_a[i].imag();                        
        }
        for(size_t i = NF; i < 2*NF; ++i) {
            wspace[i][0]=0;
            wspace[i][1]=0;
        }
        
        fftw_execute(p1); /* repeat as needed */
        for(size_t i = 0; i < 2*NF; ++i)  {
            wspace[i][0]=wspace[i][0]*wspace[i][0]+wspace[i][1]*wspace[i][1];
            wspace[i][1]=0;  
        }
        fftw_execute(p2); /* repeat as needed */
    
        for(size_t i = 0; i < NF; ++i) {
            correlated_a[i]=complex<double>(wspace[i][0],wspace[i][1])*( 1.0 / ( 2*NF * (NF -i ) ) );
        }
        
        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_free(wspace);
       
    } else {
        
        Err::Inst()->write("Correlation method not understood");
        Err::Inst()->write("scattering.correlation.method == direct, fftw");        
        throw;
        
    }
}

void AllVectorsThreadScatterDevice::infinite_correlate(std::vector<std::complex<double> >& data) {
    
    size_t NF = data.size();
          
    std::vector< std::complex<double> >& complete_a = data;
    
    complex<double> asum = 0;
    for(size_t tau = 0; tau < NF; ++tau)
    {
   		 asum += complete_a[tau];
    }

    asum /= NF;
    asum *= conj(asum);
    
    for(size_t tau = 0; tau < NF; ++tau)
    {
        complete_a[tau] = asum;
    }
}

void AllVectorsThreadScatterDevice::conjmultiply(std::vector<std::complex<double> >& data) {
    
    size_t NF = data.size();
    for(size_t n = 0; n < NF; n++)
    {
        data[n]*=conj(data[n]);
    }
}

void AllVectorsThreadScatterDevice::worker1(size_t NM,CartesianCoor3D q) {
    //size_t NMMAX = 1000;
    
    boost::threadpool::pool tp(Params::Inst()->limits.computation.threads);
    for(size_t mi = 0; mi < NM; mi+=1)
    {
//        boost::thread t(boost::bind(&AllVectorsThreadScatterDevice::scatter,this,mi) );	
        boost::threadpool::schedule(tp,boost::bind(&AllVectorsThreadScatterDevice::scatter,this,mi));
    }
       // only allow NMMAX results to reside in the internal buffer
//        while (at1.size()>=NMMAX) {
//            timer.start("sd:scatter::wait");        
//            boost::xtime xt;
//            boost::xtime_get(&xt,boost::TIME_UTC);
//            xt.nsec += 1000;
//            boost::this_thread::sleep(xt);
//            timer.stop("sd:scatter::wait");        
//        } 	    
//    }
    tp.wait();

    // push a sentinal value on the queue:
    std::pair<size_t,std::vector<complex<double> >* > sentinal(NM,NULL);
    at1.push(sentinal);

    // this feature is dead:
	//    if (Params::Inst()->scattering.center) {
    //        multiply_alignmentfactors(q);
	//    }
}


void AllVectorsThreadScatterDevice::worker2() {
    
    size_t NN = m_fqtcomm.size();
    size_t NF = p_sample->coordinate_sets.size();
    size_t NMYF = myframes.size();
    
    size_t MYRANK = m_fqtcomm.rank();
    
    // local data holders:
    std::map<size_t,size_t> qindex_recv_count;
    std::map<size_t,std::vector<std::complex<double> >* > recv_queue;
    vector<size_t> nframes(NN);
    boost::mpi::all_gather(m_fqtcomm,NMYF,&(nframes[0]));
    vector<size_t> nframe_offsets(NN);
    size_t nframe_integral=0;
    for (size_t i=0;i<NN;i++) {
        nframe_offsets[i]=nframe_integral;
        nframe_integral+=nframes[i];        
    }
    
    while ( true ) {
        
        timer.start("sd:traffic"); 
        // get a local copy of the data
        std::pair<size_t,std::vector< std::complex<double> >* > at_local_pair(0,NULL);
        std::vector< std::complex<double>  >* p_at_local = NULL;
        //if (!at1.try_pop(p_at_local)) continue;
        at1.wait_and_pop(at_local_pair);        
        
        size_t qindex = at_local_pair.first;
        p_at_local = at_local_pair.second;
        
        // test for sentinal value:
        if (p_at_local==NULL) {
            // push a sentinal value for the next worker
            at2.push(NULL);
            break;
        }

        // target node: vector_count % NN
        size_t target_node = (qindex%NN);

        std::vector<size_t> allqindexes(NN);
        boost::mpi::all_gather(m_fqtcomm,qindex,&(allqindexes[0]));
        
        std::list<boost::mpi::request> requests;
        double* p_at_frames = (double*) &((*p_at_local)[0]);
        requests.push_back(m_fqtcomm.isend(target_node,0,p_at_frames,2*NMYF));
 
        // loop through all qindexes and post a recv for each qindex which matches Myrank as a target
        for (size_t i=0;i<NN;i++) {
            size_t this_qindex = allqindexes[i];
            
            size_t this_target_node = (this_qindex%NN);
            if (MYRANK==this_target_node) {
                std::vector<std::complex<double> >* p_at_allframes = NULL;
                
                // reserve data on local queue if qindex is new
                if (recv_queue.find(this_qindex)==recv_queue.end()) {
                    p_at_allframes = new std::vector<std::complex<double> >(NF);
                    recv_queue[this_qindex]=p_at_allframes; // push the pointer
                } else {
                    p_at_allframes = recv_queue[this_qindex];
                }
                
                // select the region the data is recv
                double* pp_at_allframes = (double*) &((*p_at_allframes)[nframe_offsets[i]]);
                // post a blocking recv                
                m_fqtcomm.recv(i,0,pp_at_allframes,2*nframes[i]);                
                
                qindex_recv_count[this_qindex]++;
            }
        }
        
        
        boost::mpi::wait_all(requests.begin(),requests.end());
	    // need to free p_at_local, since we gained ownership
        delete p_at_local;
        
        std::queue<size_t> qindex_finished;
        // now scan for completely recv data
        for (std::map<size_t,size_t>::iterator qiter=qindex_recv_count.begin();qiter!=qindex_recv_count.end();qiter++) {
            if (qiter->second==NN) {
                at2.push(recv_queue[qiter->first]);
                qindex_finished.push(qiter->first);
            }            
        }
        while (!qindex_finished.empty()) {
            size_t qindex = qindex_finished.front();
            recv_queue.erase(qindex);
            qindex_recv_count.erase(qindex);            
            qindex_finished.pop();
        }
	}
}

void AllVectorsThreadScatterDevice::worker3() {
    
    //size_t NF = p_sample->coordinate_sets.size();
    while ( true ) {
        // get a local copy of the data
        std::vector< std::complex<double> >* p_at_local=NULL;
        //if (!at2.try_pop(p_at_local)) continue;
        at2.wait_and_pop(p_at_local); 
              
        if (p_at_local==NULL) {
             break;
        }
        
        timer.start("sd:correlation"); 

		// correlate or sum up
	    if (Params::Inst()->scattering.correlation.type=="time") {

    	    timer.start("sd:correlate");	                    
            correlate(*p_at_local);
    	    timer.stop("sd:correlate");	        

    	} else if (Params::Inst()->scattering.correlation.type=="infinite-time") {

            timer.start("sd:correlate");	                    
            infinite_correlate(*p_at_local);
            timer.stop("sd:correlate");	        
        
		} else {

    	    timer.start("sd:conjmultiply");	                                
            conjmultiply(*p_at_local);
    	    timer.stop("sd:conjmultiply");	                    
		}
		
        for(size_t n=0;n<at3.size();n++) {
            at3[n]+=(*p_at_local)[n];
        }
        timer.stop("sd:correlation"); 

        delete p_at_local;
	}	
}


void AllVectorsThreadScatterDevice::next() {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	m_current_qvector+=1;    
	
    ofstream monitorfile("progress.data",ios::binary);
    monitorfile.seekp(m_scattercomm.rank()*sizeof(float));
        float progress = m_current_qvector*1.0/m_qvectorindexpairs.size();
    monitorfile.write((char*) &progress,sizeof(float)); // initialize everything to 0
    monitorfile.close();
}

double AllVectorsThreadScatterDevice::progress() {
    return m_current_qvector*1.0/m_qvectorindexpairs.size();    
}

size_t AllVectorsThreadScatterDevice::status() {
    return m_current_qvector/m_qvectorindexpairs.size();
}

void AllVectorsThreadScatterDevice::write() {
 // do a scatter_comm operation to negotiate write outs
 for(int i = 0; i < m_scattercomm.size(); ++i)
 {
     if (m_scattercomm.rank()==i) {
    	if ((m_fqtcomm.rank()==0) && m_writeflag) {
             H5FQTInterface::store(m_fqt_filename,m_qvectorindexpairs[m_current_qvector].first,m_spectrum);
             m_writeflag = false;
    	}
    }
    m_scattercomm.barrier();
 }
}


void AllVectorsThreadScatterDevice::compute() {

	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	CartesianCoor3D q=m_qvectorindexpairs[m_current_qvector].second;

    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	// we need 3 worker threads:
	// one for the calculating of the complex amplitudes
	// one for the network communication / transposition of the complex amplitudes
	// one for the calculation of the correlation

    // start worker1
    // start worker2
    // start worker3

    
	// blocking factor: max(nn,nq)
    size_t NN = m_fqtcomm.size();
    size_t NM = get_numberofmoments();
    size_t NF = p_sample->coordinate_sets.size();

    m_spectrum.clear();
    at3.clear();
    m_spectrum.resize(NF,0);
    at3.resize(NF,0);
    
    boost::thread t1(boost::bind(&AllVectorsThreadScatterDevice::worker1, this, NM,q));
    boost::thread t2(boost::bind(&AllVectorsThreadScatterDevice::worker2, this));    
    boost::thread t3(boost::bind(&AllVectorsThreadScatterDevice::worker3, this));

    // t3 is the last to go
    t3.join();

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

    m_writeflag=true;
}


void AllVectorsThreadScatterDevice::scatter(size_t moffset) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;
   
   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = myframes.size();

   std::vector<complex<double> >* p_a = new std::vector<complex<double> >(NMYF,0);
   
   string target = Params::Inst()->scattering.target;
   size_t NOA = p_sample->atoms.selections[target].indexes.size();

   CartesianCoor3D& q = qvectors[moffset];

   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       boost::mutex::scoped_lock lock(scatter_mutex);
       timer.start("sd:fs:f:ld");	
       CoordinateSet& cs = *csets[fi]; 
       timer.stop("sd:fs:f:ld");	
       lock.unlock();
               
       double Ar =0; 
       double Ai = 0;
       double x,y,z;
       double qx = q.x;
       double qy = q.y;
       double qz = q.z;
       double esf,p;
       
	   for(size_t j = 0; j < NOA; ++j) {   		
	       
	       esf = sfs[j];
           x = cs.c1[j];
   	       y = cs.c2[j];
   	       z = cs.c3[j];
       
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

size_t AllVectorsThreadScatterDevice::get_numberofmoments() {
    return qvectors.size();
}
