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
#include "scatter_devices/self_vectors.hpp"

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
#include "scatter_devices/particle_trajectory.hpp"


using namespace std;

SelfVectorsScatterDevice::SelfVectorsScatterDevice(boost::mpi::communicator& thisworld, Sample& sample){

	p_thisworldcomm = &thisworld;
	p_sample = &sample; // keep a reference to the sample

	string target = Params::Inst()->scattering.target;
	
	size_t NN = thisworld.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	size_t rank = thisworld.rank();

	EvenDecompose edecomp(NF,NN);

	if (NN>NA) {
		Err::Inst()->write("Not enough atoms for decomposition");
		throw;
	}

	size_t minatoms = NA/NN;	
	size_t leftoveratoms = NA % NN;	

//	CoordinateSet* p_cs =NULL;

	size_t* p_indexes; size_t* p_myindexes;
	// construct particle_trajectories first
	vector<size_t> selectionindexes(NA);
	if (rank==0) {
//		sample.frames.load(0,sample.atoms);
//		p_cs = new CoordinateSet(sample.frames.current(),sample.atoms.selections[target]);
		// don't use the real atom indexes, but the selection indexes instead, this will work well with the scatterfactors class			
//		p_indexes = &(sample.atoms.selections[target][0]); 
		for(size_t i = 0; i < selectionindexes.size(); ++i) selectionindexes[i]=i;
		p_indexes = &(selectionindexes[0]);
	}

	vector<size_t> myindexes(minatoms); p_myindexes = &(myindexes[0]);
	
	timer.start("sd:assign");	
	// this is an atom decomposition!
	boost::mpi::scatter(thisworld,p_indexes,p_myindexes,minatoms,0);	
	for(size_t i = 0; i < minatoms; ++i)
	{
		particle_trajectories.push_back(ParticleTrajectory(myindexes[i]) );
	}

	if (leftoveratoms>0) {
		
		if (rank==0) {
			p_indexes = &(sample.atoms.selections[target].indexes[minatoms]);
		}

		size_t lastindex; p_myindexes = &(lastindex);
		boost::mpi::scatter(thisworld,p_indexes,p_myindexes,1,0);	
		if (rank<leftoveratoms) {
			particle_trajectories.push_back(ParticleTrajectory(lastindex) );		
		}	
	}
	timer.stop("sd:assign");	

	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);

	if (Params::Inst()->scattering.center) {
    	p_sample->coordinate_sets.add_postalignment(target,"center");		    
	}

	timer.start("sd:data");	
	for(size_t r = 0; r < NN; ++r)
	{

		vector<size_t> assigned_frames = edecomp.indexes_for(r);
		
		for(size_t i = 0; i < assigned_frames.size(); ++i)
		{
			double  *p_xc,*p_yc,*p_zc;	
			
			CoordinateSet* p_cs =NULL;
			if (rank==r) {
                p_cs = &( p_sample->coordinate_sets.load(assigned_frames[i]) );
                                
				p_xc = &(p_cs->c1[0]);
				p_yc = &(p_cs->c2[0]);
				p_zc = &(p_cs->c3[0]);
							
			}
			
			vector<double> myxc(minatoms);
			vector<double> myyc(minatoms);
			vector<double> myzc(minatoms);

			boost::mpi::scatter(thisworld,p_xc,&(myxc[0]),minatoms,r);
			boost::mpi::scatter(thisworld,p_yc,&(myyc[0]),minatoms,r);
			boost::mpi::scatter(thisworld,p_zc,&(myzc[0]),minatoms,r);	

			for(size_t pi = 0; pi < minatoms; ++pi)
			{
				particle_trajectories[pi].push_back(CartesianCoor3D(myxc[pi],myyc[pi],myzc[pi]));
			}		

			if (rank==r) {
				p_xc = &(p_cs->c1[minatoms]);
				p_yc = &(p_cs->c2[minatoms]);
				p_zc = &(p_cs->c3[minatoms]);
			}

	        // we have to exchange alignment information
            // b/c the individual nodes don't have all the alignment vectors at once
            // but they are needed
	        if (Params::Inst()->scattering.center) {
                vector<CartesianCoor3D> myavectors;
                if (rank==r) {
                	myavectors = p_sample->coordinate_sets.get_postalignmentvectors(assigned_frames[i]);
			    }
    		    boost::mpi::broadcast(thisworld,myavectors,r);
                m_all_postalignmentvectors[assigned_frames[i]]=myavectors;
	        }
	        
			double mylastx; 
			double mylasty; 
			double mylastz; 
			
			boost::mpi::scatter(thisworld,p_xc,&mylastx,1,r);
			boost::mpi::scatter(thisworld,p_yc,&mylasty,1,r);
			boost::mpi::scatter(thisworld,p_zc,&mylastz,1,r);	

			if (rank<leftoveratoms) {
				particle_trajectories[minatoms].push_back(CartesianCoor3D(mylastx,mylasty,mylastz));
			}
				
		}
	}
	timer.stop("sd:data");	
	
	p_asingle = new vector< complex<double> >; // initialize with no siz
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);
	
}

void SelfVectorsScatterDevice::scatter(size_t ai, size_t mi) {

	ParticleTrajectory& thisparticle = particle_trajectories[ai];
	// this is broken <-- revise this!!!
	double s = scatterfactors.get(thisparticle.atom_selection_index());
	
    size_t NF = p_sample->coordinate_sets.size();
	
    p_asingle->resize(NF,0);	

    double qx = qvectors[mi].x;
    double qy = qvectors[mi].y;
    double qz = qvectors[mi].z;
    
	for(size_t j = 0; j < NF; ++j)
	{
		double x1 = thisparticle.x[j];
		double y1 = thisparticle.y[j];
		double z1 = thisparticle.z[j];	
		double p1 = x1*qx+y1*qy+z1*qz;
        
        double sp1 = sin(p1);
        double cp1 = cos(p1);
        (*p_asingle)[j] = s*complex<double>(cp1,sp1);
	}
	
}

void SelfVectorsScatterDevice::multiply_alignmentfactors(size_t mi) {

    CartesianCoor3D& q = qvectors[mi];
    
    for(size_t j = 0; j < p_asingle->size(); ++j)
    {
        // a2 index is frame number!
        vector<CartesianCoor3D>& avectors = m_all_postalignmentvectors[j];
        CartesianCoor3D& bigR = *(avectors.rbegin());
        complex<double> factor = exp(complex<double>(0,qvectors[mi]*bigR));	        

        (*p_asingle)[j] = (*p_asingle)[j] * factor;
    }    
}


void SelfVectorsScatterDevice::correlate() {
    if (p_asingle->size()<1) return;
    
    size_t NF = p_sample->coordinate_sets.size();
    
    std::vector< std::complex<double> >* p_correlated_a = new std::vector< std::complex<double> >;
    p_correlated_a->resize(NF);
      
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

void SelfVectorsScatterDevice::execute(CartesianCoor3D q) {
	
    init(q);

	timer.start("sd:sf:update");
	scatterfactors.update(q); // scatter factors only dependent on length of q, hence we can do it once before the loop
	timer.stop("sd:sf:update");
	
	// blocking factor: max(nn,nq)
    long NN = p_thisworldcomm->size();
    long NM = get_numberofmoments();
    long NF = p_sample->coordinate_sets.size();
    
	m_spectrum.resize(NF,0);

    timer.start("sd:compute");
    // correlate or sum up
    if (Params::Inst()->scattering.correlation.type=="time") {

        for (long ai = 0; ai < particle_trajectories.size(); ai++) {
            for(long mi = 0; mi < NM; mi++)
            {        
                // compute a block of q vectors:
    	        timer.start("sd:scatterblock");
    	        scatter(ai,mi);	// fills p_asingle
    	        timer.stop("sd:scatterblock");            

    	        // operate on (*p_asingle)
    	        if (Params::Inst()->scattering.center) {
                    multiply_alignmentfactors(mi);
    	        }		

    	        timer.start("sd:correlate");	                    
                correlate();
    	        timer.stop("sd:correlate");	        

                vector<complex<double> >& spectrum = (*p_asingle);
        	    for(size_t j = 0; j < spectrum.size(); ++j)
        	    {
        	    	m_spectrum[j] += spectrum[j];
        	    }	
        	}				
    	}
    	
    } else {
        // if not time correlated, the conjmultiply negates phase information
        // this simplifies formulas
 
        p_asingle->resize(NF,0);
 
        for (long ai = 0; ai < particle_trajectories.size(); ai++) {
            ParticleTrajectory& thisparticle = particle_trajectories[ai];
        	// this is broken <-- revise this!!!
        	double s = scatterfactors.get(thisparticle.atom_selection_index());
            for(size_t fi = 0; fi < NF; ++fi)
            {
                m_spectrum[fi] += s*s;
            }
        }
    	
    }
    timer.start("sd:compute");    
    	
	
	// these functions operate on m_spectrum:
	
    timer.start("sd:norm");	                                		
    norm();
    timer.stop("sd:norm");	                    

	// finally assemble results on the head node:
    timer.start("sd:gather_sum");	                    
    gather_sum(); // head node has everything in a (vector< complex<double> >)
    timer.stop("sd:gather_sum");	        

}


void SelfVectorsScatterDevice::gather_sum() {
    size_t NN = p_thisworldcomm->size();
    size_t NF = p_sample->coordinate_sets.size();
        
    if (NF<1) return;
    
    if (p_thisworldcomm->rank()==0) {
        vector< complex<double> > received_a(NF);
        double* p_a_double = (double*) &(received_a.at(0));   
        vector< complex<double> >& a = m_spectrum;
        
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
        double* p_a_double = (double*) &(m_spectrum[0]);   
        p_thisworldcomm->send(0,0,p_a_double,2*NF);
        m_spectrum.clear();
    }

}


vector<complex<double> >& SelfVectorsScatterDevice::get_spectrum() {
	return m_spectrum;
}


void SelfVectorsScatterDevice::init(CartesianCoor3D& q) {
    
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

size_t SelfVectorsScatterDevice::get_numberofmoments() {
    return qvectors.size();
}

void SelfVectorsScatterDevice::norm() {
    size_t NV = p_asingle->size();
    for(size_t i = 0; i < NV; ++i)
    {
        (*p_asingle)[i] /= qvectors.size();
    }
}

// end of file