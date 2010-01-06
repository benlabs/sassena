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
#include <fftw3.h>
#include "control.hpp"
#include "sample.hpp"

using namespace std;

AllScatterDevice::AllScatterDevice(boost::mpi::communicator& thisworld, Sample& sample) {

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;

	size_t NN = thisworld.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	size_t rank = thisworld.rank();

	EvenDecompose edecomp(NF,NN);
	
	myframes = edecomp.indexes_for(rank);
    // blocking factor:
    // max(nn,nq)

	if (thisworld.rank()==0) {	

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
			Warn::Inst()->write(string("Insufficient Buffer size for coordinate sets. This is a HUGE bottleneck for performance!"));
			Warn::Inst()->write(string("Need at least: ")+to_s(memusage_allcs)+string(" bytes"));
			Warn::Inst()->write(string("Configuration Parameter: limits.memory.coordinate_sets"));
		}		
	}
	
    p_a = new vector< vector< complex<double> > >; // initialize with no size
    p_asingle = new vector< complex<double> >; // initialize with no size

	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

    // overwrite representation style in subclass !!
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	if (Params::Inst()->scattering.center) {
    	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	}

	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
}

AllScatterDevice::~AllScatterDevice() {
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
    size_t NN = p_thisworldcomm->size();
    // first gather number of frames per process:
    vector<size_t> nframes(NN);
    boost::mpi::all_gather(*p_thisworldcomm,NMYF,&(nframes[0]));

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
    if (p_thisworldcomm->rank()<blocknq) {
        for(size_t i = 0; i < NN; ++i)
        {
            double* p_newa_frames = (double*) &(newa[i][0]);
            requests.push_back( p_thisworldcomm->irecv(i,0,p_newa_frames,2*nframes[i]) );
        }        
    }

    for(size_t i = 0; i < blocknq; ++i)
    {
        double* p_a_frames = (double*) &((*p_a)[i][0]);        
        requests.push_back( p_thisworldcomm->isend(i,0,p_a_frames,2* (*p_a)[i].size() ) );
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
    size_t NN = p_thisworldcomm->size();
    size_t NF = p_sample->coordinate_sets.size();
        
    if (NF<1) return;
    
    if (p_thisworldcomm->rank()==0) {
        vector< complex<double> > received_a(NF);
        double* p_a_double = (double*) &(received_a.at(0));   
        vector< complex<double> >& a = (*p_asingle);
        
        // only receive from other nodes
        for(size_t i = 1; i < NN; ++i)
        {
            p_thisworldcomm->recv(i,0,p_a_double,2*NF);

            for(size_t i = 0; i < NF; ++i)
            {
                a[i] += received_a[i];
            }
        }

    } else {
        double* p_a_double = (double*) &((*p_asingle)[0]);   
        p_thisworldcomm->send(0,0,p_a_double,2*NF);
        p_asingle->clear();
    }

}

void AllScatterDevice::gather_cat() {
    size_t NN = p_thisworldcomm->size();
    size_t NF = p_sample->coordinate_sets.size();
            
    if (NF<1) return;

    // send frames of j'th column to j'th node
    size_t NMYF = myframes.size();
    // first gather number of frames per process:
    vector<size_t> nframes(NN);
    boost::mpi::all_gather(*p_thisworldcomm,NMYF,&(nframes[0]));
    
    if (p_thisworldcomm->rank()==0) {
        vector< complex<double> >& a = (*p_asingle);
        
        // only receive from other nodes
        for(size_t i = 1; i < NN; ++i)
        {
            vector< complex<double> > received_a(nframes[i]);
            double* p_a_double = (double*) &(received_a[0]);   
            p_thisworldcomm->recv(i,0,p_a_double,2*nframes[i]);
            
            size_t thisNF = nframes[i];
            for(size_t j = 0; j < thisNF; ++j)
            {
                a.push_back( received_a[j] );
            }
        }

    } else {
        double* p_a_double = (double*) &(p_asingle->at(0));   
        p_thisworldcomm->send(0,0,p_a_double,2*NMYF);
        p_asingle->clear();
    }
    
}

void AllScatterDevice::execute(CartesianCoor3D q) {

    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	// blocking factor: max(nn,nq)
    long NN = p_thisworldcomm->size();
    long NM = get_numberofmoments();
    long NF = p_sample->coordinate_sets.size();
    long NMYF = myframes.size();

    long NMBLOCK = (NN<NM) ? NN : NM;
    
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
}

vector<complex<double> >& AllScatterDevice::get_spectrum() {
	return m_spectrum;
}
