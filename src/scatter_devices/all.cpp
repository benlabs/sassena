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
#include "scatter_devices/all.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

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

AllScatterDevice::AllScatterDevice(	boost::mpi::communicator all_comm,
		boost::mpi::communicator partition_comm,
		Sample& sample,
		vector<pair<size_t,CartesianCoor3D> > QIV,
		string fqt_filename)
{
	m_qvectorindexpairs= QIV;
	m_current_qvector =0;
	m_fqt_filename = fqt_filename;

	m_scattercomm = all_comm;
	m_fqtcomm = partition_comm;
    m_writeflag = false;

	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t NN = m_fqtcomm.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	size_t rank = m_fqtcomm.rank();
	
	//////////////////////////////////////////////////
	// Assignment of frames to compute
	//////////////////////////////////////////////////
	EvenDecompose edecomp(NF,NN);
	myframes = edecomp.indexes_for(rank);
	
	
    // blocking factor:
    // max(nn,nq)
	
    p_a = new vector< vector< complex<double> > >; // initialize with no size
    p_asingle = new vector< complex<double> >; // initialize with no size

	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);
    // overwrite representation style in subclass !!
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	if (Params::Inst()->scattering.center) {
    	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	}
	
    for (size_t mf=0;mf<myframes.size();mf++) {
        csets.push_back(p_sample->coordinate_sets.load(myframes[mf]));
    }
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
}

AllScatterDevice::~AllScatterDevice() {
    for (size_t ci=0;csets.size();ci++) {
        delete csets[ci];
    }
    
    delete p_a;
    delete p_asingle;
}

// scatter_frame implemented by concrete class
// scatter_frames implemented by concrete class

// init_averaging implemented by concrete class
// get_numberofmoments implemented by concrete class

void AllScatterDevice::multiply_alignmentfactors(CartesianCoor3D q) {
    
    for(size_t i = 0; i < p_a->size(); ++i)
	{
        vector< complex<double> >& frames = p_a->at(i);
        if (frames.size()<1) continue;
        
        size_t iframe = myframes[i];
        vector<CartesianCoor3D> avectors = p_sample->coordinate_sets.get_postalignmentvectors(iframe);
        CartesianCoor3D& bigR = *(avectors.rbegin());
        complex<double> factor = exp(complex<double>(0,q*bigR));
        
        size_t NF = p_a->at(i).size();
	    for(size_t j = 0; j < NF; ++j)
	    {
            frames[j] = frames[j] * factor;
	    }
    }
    
}

// transpose replaces the scattering matrix a (firstq,secondframes) with a transposed version (firstframes,secondq)
void AllScatterDevice::exchange() {
    
    // send frames of j'th column to j'th node
    size_t blocknq = p_a->size();
    if (blocknq<1) return;
    
    size_t NMYF = myframes.size();
    size_t NN = m_fqtcomm.size();
    size_t FQT_RANK = m_fqtcomm.rank();
    // first gather number of frames per process:
    vector<size_t> nframes(NN);
    boost::mpi::all_gather(m_fqtcomm,NMYF,&(nframes[0]));

    // this vector is indexed by the process
    std::vector< std::vector< std::complex<double> > > newa(NN);
    for(size_t i = 0; i < NN; ++i)
    {
        newa[i].resize(nframes[i]);
    } 
    
    // this will be scoped to the relevant q vector
    std::list< boost::mpi::request > requests;

    // exchange real part 
    // pre-post receives
    if (FQT_RANK<blocknq) {
        for(size_t i = 0; i < NN; ++i)
        {
            double* p_newa_frames = (double*) &(newa[i][0]);
            requests.push_back( m_fqtcomm.irecv(i,0,p_newa_frames,2*nframes[i]) );
        }        
    }

    for(size_t i = 0; i < blocknq; ++i)
    {
        double* p_a_frames = (double*) &((*p_a)[i][0]);        
        requests.push_back( m_fqtcomm.isend(i,0,p_a_frames,2* (*p_a)[i].size() ) );
    }

    boost::mpi::wait_all(requests.begin(),requests.end());
    
    // reconstruct the time sequence of scattering amplitudes
    // process i has i'th qvector
    
    std::vector< std::complex<double> >* p_newcompletea = new std::vector< std::complex<double> >;
    
    vector< complex<double> >& frames = *p_newcompletea;
    for(size_t i = 0; i < NN; ++i)
    {
        for(size_t j = 0; j < nframes[i]; ++j)
        {
            frames.push_back(newa[i][j]);
        }
    }
    
    // switch a
    p_a->clear();
    delete p_asingle;
    p_asingle = p_newcompletea;
}

void AllScatterDevice::correlate() {
    if (p_asingle->size()<1) return;
    
    size_t NF = p_sample->coordinate_sets.size();
    
    std::vector< std::complex<double> >* p_correlated_a = new std::vector< std::complex<double> >;
    p_correlated_a->assign(NF,0);
      
    std::vector< std::complex<double> >& complete_a = (*p_asingle);
    std::vector< std::complex<double> >& correlated_a = (*p_correlated_a);
    
    if (Params::Inst()->scattering.correlation.method=="direct") {
        
        // direct
        for(size_t tau = 0; tau < NF; ++tau)
        {
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
        
    // change meaning of asingle
    delete p_asingle;
    p_asingle = p_correlated_a;
}


void AllScatterDevice::infinite_correlate() {
    if (p_asingle->size()<1) return;
    
    size_t NF = p_sample->coordinate_sets.size();
          
    std::vector< std::complex<double> >& complete_a = (*p_asingle);
    
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

void AllScatterDevice::conjmultiply() {
    
    size_t NQBLOCK = p_a->size();
    for(size_t qi = 0; qi < NQBLOCK; ++qi)
    {
        vector<std::complex<double> >& values = (*p_a)[qi];
        size_t NV = values.size();
        for(size_t vi = 0; vi < NV; ++vi)
        {
            values[vi]*=conj(values[vi]);
        }
    }
}


void AllScatterDevice::average_correlate() {
    if (p_asingle->size()<1) return;
    
    size_t NF = p_sample->coordinate_sets.size();
          
    std::vector< std::complex<double> >& complete_a = (*p_asingle);
    
    complex<double> asum = 0;
    for(size_t tau = 0; tau < NF; ++tau)
    {
   		 asum += complete_a[tau];
    }

    asum /= NF;

    for(size_t tau = 0; tau < NF; ++tau)
    {
        complete_a[tau] = asum;
    }
}

void AllScatterDevice::sum() {
    
    size_t NMYF = myframes.size();

    p_asingle->assign(NMYF,0);

    vector< complex<double> >& newa = (*p_asingle);
    for(size_t qi = 0; qi < p_a->size(); ++qi)
    {
        vector< complex<double> >& a = (*p_a)[qi];
        for(size_t fi = 0; fi < NMYF; ++fi)
        {
            newa[fi]+=a[fi];
        }
    }
    
    // free space
    p_a->clear();
}

void AllScatterDevice::gather_sum() {
    size_t NN = m_fqtcomm.size();
    size_t NF = p_sample->coordinate_sets.size();
        
    if (NF<1) return;
    
    if (m_fqtcomm.rank()==0) {
        vector< complex<double> > received_a(NF);
        double* p_a_double = (double*) &(received_a.at(0));   
        vector< complex<double> >& a = (*p_asingle);
        
        // only receive from other nodes
        for(size_t i = 1; i < NN; ++i)
        {
            m_fqtcomm.recv(i,0,p_a_double,2*NF);

            for(size_t i = 0; i < NF; ++i)
            {
                a[i] += received_a[i];
            }
        }

    } else {
        double* p_a_double = (double*) &((*p_asingle)[0]);   
        m_fqtcomm.send(0,0,p_a_double,2*NF);
        p_asingle->clear();
    }

}

void AllScatterDevice::gather_cat() {
    size_t NN = m_fqtcomm.size();
    size_t NF = p_sample->coordinate_sets.size();
            
    if (NF<1) return;

    // send frames of j'th column to j'th node
    size_t NMYF = myframes.size();
    // first gather number of frames per process:
    vector<size_t> nframes(NN);
    boost::mpi::all_gather(m_fqtcomm,NMYF,&(nframes[0]));
    
    if (m_fqtcomm.rank()==0) {
        vector< complex<double> >& a = (*p_asingle);
        
        // only receive from other nodes
        for(size_t i = 1; i < NN; ++i)
        {
            vector< complex<double> > received_a(nframes[i]);
            double* p_a_double = (double*) &(received_a[0]);   
            m_fqtcomm.recv(i,0,p_a_double,2*nframes[i]);
            
            size_t thisNF = nframes[i];
            for(size_t j = 0; j < thisNF; ++j)
            {
                a.push_back( received_a[j] );
            }
        }

    } else {
        double* p_a_double = (double*) &(p_asingle->at(0));   
        m_fqtcomm.send(0,0,p_a_double,2*NMYF);
        p_asingle->clear();
    }
}


void AllScatterDevice::next() {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	m_current_qvector+=1;    
	
    ofstream monitorfile("progress.data",ios::binary);
    monitorfile.seekp(m_scattercomm.rank()*sizeof(float));
    float progress = m_current_qvector*1.0/m_qvectorindexpairs.size();
    monitorfile.write((char*) &progress,sizeof(float)); // initialize everything to 0
    monitorfile.close();
}

double AllScatterDevice::progress() {
    return m_current_qvector*1.0/m_qvectorindexpairs.size();    
}

size_t AllScatterDevice::status() {
    return m_current_qvector/m_qvectorindexpairs.size();
}


void AllScatterDevice::compute() {

	if (m_current_qvector>=m_qvectorindexpairs.size()) return;
	CartesianCoor3D q=m_qvectorindexpairs[m_current_qvector].second;

    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	// blocking factor: max(nn,nq)
    long NN = m_fqtcomm.size();
    long NM = get_numberofmoments();
    long NF = p_sample->coordinate_sets.size();
    //long NMYF = myframes.size();

    long NMBLOCK = (NN<NM) ? NN : NM;
//    NMBLOCK = 1;
    
    m_spectrum.assign(NF,0);
    for(long mi = 0; mi < NM; mi+=NMBLOCK)
    {        

        // compute a block of q vectors:
	    timer.start("sd:scatterblock");
	    scatter(mi,std::min(NMBLOCK,NM-mi));	
	    timer.stop("sd:scatterblock");
	    if (Params::Inst()->scattering.center) {
            multiply_alignmentfactors(q);
	    }		
                
		// correlate or sum up
	    if (Params::Inst()->scattering.correlation.type=="time") {

	        // transpose data
    	    timer.start("sd:exchange");	        
            exchange();
    	    timer.stop("sd:exchange");	        

    	    timer.start("sd:correlate");	                    
            correlate();
    	    timer.stop("sd:correlate");	        

    	    timer.start("sd:gather_sum");	                    
            gather_sum(); // head node has everything in a (vector< complex<double> >)
    	    timer.stop("sd:gather_sum");
    	    	        
    	} else if (Params::Inst()->scattering.correlation.type=="infinite-time") {

            // transpose data
            timer.start("sd:exchange");	        
            exchange();
            timer.stop("sd:exchange");	        
        
            timer.start("sd:correlate");	                    
            infinite_correlate();
            timer.stop("sd:correlate");	        
        
            timer.start("sd:gather_sum");	                    
            gather_sum(); // head node has everything in a (vector< complex<double> >)
            timer.stop("sd:gather_sum");	        
        
		} else {

    	    timer.start("sd:conjmultiply");	                                
            conjmultiply();
    	    timer.stop("sd:conjmultiply");	                    

    	    timer.start("sd:sum");	                                
            sum();
    	    timer.stop("sd:sum");	                    

    	    timer.start("sd:gather_cat");	                                
            gather_cat();
    	    timer.stop("sd:gather_cat");	                    

		}

	    timer.start("sd:norm");	                                		
        norm();
	    timer.stop("sd:norm");	                    
					
        vector<complex<double> >& spectrum = (*p_asingle);
    	for(size_t j = 0; j < spectrum.size(); ++j)
    	{
    		m_spectrum[j] += spectrum[j];
    	}	
    						
	}

    m_writeflag=true;
}

void AllScatterDevice::write() {
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
