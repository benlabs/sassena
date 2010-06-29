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

// other headers
#include "math/coor3d.hpp"
#include "math/smath.hpp"
#include "decomposition/decompose.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;

AllVectorsThreadScatterDevice::AllVectorsThreadScatterDevice(boost::mpi::communicator& thisworld, Sample& sample) {

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
	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

    // overwrite representation style in subclass !!
	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	
	//if (Params::Inst()->scattering.center) {
    //	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	//}

	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
	
    size_t NMYF = myframes.size();
	
    size_t NM = 1;
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        NM = Params::Inst()->scattering.average.orientation.vectors.size();
    }
    size_t NMBLOCK = (NN<NM) ? NN : NM;
    
	if (thisworld.rank()==0) {	

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
    size_t NMMAX = 1000;
    for(size_t mi = 0; mi < NM; mi+=1)
    {
        timer.start("sd:scatter");
        scatter(mi);	
        timer.stop("sd:scatter");        
        
       // only allow NMMAX results to reside in the internal buffer
        while (at1.size()>=NMMAX) {
            timer.start("sd:scatter::wait");        
            boost::xtime xt;
            boost::xtime_get(&xt,boost::TIME_UTC);
            xt.nsec += 1000;
            boost::this_thread::sleep(xt);
            timer.stop("sd:scatter::wait");        
        } 	    
    }
    
    worker1_done_flag = true;
    // this feature is dead:
	//    if (Params::Inst()->scattering.center) {
    //        multiply_alignmentfactors(q);
	//    }
}

void AllVectorsThreadScatterDevice::worker2() {
    
    size_t NN = p_thisworldcomm->size();
    size_t NF = p_sample->coordinate_sets.size();
    size_t NMYF = myframes.size();
    
    size_t MYRANK = p_thisworldcomm->rank();
    size_t vector_count = 0;
    
    while ( !(worker1_done_flag && at1.empty()) ) {
        
        timer.start("sd:traffic"); 
        // get a local copy of the data
        std::vector< std::complex<double> >* p_at_local;
        //if (!at1.try_pop(p_at_local)) continue;
        at1.wait_and_pop(p_at_local);


        // target node: vector_count % NN
        size_t target_node = (vector_count%NN);
        vector_count++;	

        vector<size_t> nframes(NN);
        boost::mpi::all_gather(*p_thisworldcomm,NMYF,&(nframes[0]));
        
        if (MYRANK==target_node) {
            
            std::vector<std::complex<double> >* p_at_allframes = new std::vector<std::complex<double> >(NF);
            
            // receive all values
            size_t ntotal=0;
            for (size_t n=0;n<NN;n++) {
                double* pp_at_allframes = (double*) &((*p_at_allframes)[ntotal]);
                if (MYRANK!=n) {
                    p_thisworldcomm->recv(n,0,pp_at_allframes,2*nframes[n]);
                } else {
                    for(size_t nn=0;nn<nframes[n];nn++) {
                        (*p_at_allframes)[ntotal+nn]=(*p_at_local)[nn];                        
                    }
                }
                ntotal+=nframes[n];
            }
            
            // push them off 
            at2.push(p_at_allframes);
                            
        } else {
            double* p_at_frames = (double*) p_at_local;
            p_thisworldcomm->send(target_node,0,p_at_frames,2*NMYF);                
        }
        timer.stop("sd:traffic"); 

	    // need to free p_at_local, since we gained ownership
        delete p_at_local;
        
	}
	
    worker2_done_flag=true;
}

void AllVectorsThreadScatterDevice::worker3() {
    
    size_t NF = p_sample->coordinate_sets.size();

    while ( !(worker2_done_flag && at2.empty()) ) {
            
        // get a local copy of the data
        std::vector< std::complex<double> >* p_at_local;
        //if (!at2.try_pop(p_at_local)) continue;
        at2.wait_and_pop(p_at_local);

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
	
    worker3_done_flag=true;
}

void AllVectorsThreadScatterDevice::execute(CartesianCoor3D q) {

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
    long NN = p_thisworldcomm->size();
    long NM = get_numberofmoments();
    long NF = p_sample->coordinate_sets.size();

    m_spectrum.assign(NF,0);
    at3.resize(NF,0);
    
    worker1_done_flag = false;
    worker2_done_flag = false;
    worker3_done_flag = false;

    boost::thread t1(boost::bind(&AllVectorsThreadScatterDevice::worker1, this, NM,q));
    
    boost::thread t2(boost::bind(&AllVectorsThreadScatterDevice::worker2, this));
    boost::thread t3(boost::bind(&AllVectorsThreadScatterDevice::worker3, this));
	    
	// wait for t4 to finish
    t3.join();

    vector<complex<double> > at_local(NF,0);
    
    double* p_at3 = (double*) &(at3[0]);
    double* p_atlocal = (double*) &(at_local[0]);
    
    
    if (NN>1) {
        boost::mpi::reduce(*p_thisworldcomm,p_at3,2*NF,p_atlocal,std::plus<double>() ,0);        
        for(size_t j = 0; j < at_local.size(); ++j)
        {
            m_spectrum[j] += std::complex<double>((at_local[j].real())/NM,(at_local[j].imag())/NM);
        }         			
    } else {
        for(size_t j = 0; j < at3.size(); ++j)
        {
            m_spectrum[j] += std::complex<double>((at3[j].real())/NM,(at3[j].imag())/NM);
        }         			
    }
   
}

vector<complex<double> >& AllVectorsThreadScatterDevice::get_spectrum() {
	return m_spectrum;
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
       timer.start("sd:fs:f:ld");	
       CoordinateSet& cs = p_sample->coordinate_sets.load(myframes[fi]); 
       timer.stop("sd:fs:f:ld");	
               
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
	}
	
    at1.push(p_a);
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
			double qparallel_l = qparallel.length();

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
