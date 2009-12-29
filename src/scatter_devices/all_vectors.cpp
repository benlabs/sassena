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
#include "scatter_devices/all_vectors.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

// other headers
#include "analysis.hpp"
#include "coor3d.hpp"
#include "decompose.hpp"
#include "fftw/fftw++.h"
#include <fftw3.h>
#include "control.hpp"
#include "sample.hpp"
#include "smath.hpp"
#include "particle_trajectory.hpp"
#include "vector_unfold.hpp"

using namespace std;


AllVectorsScatterDevice::AllVectorsScatterDevice(boost::mpi::communicator& thisworld, Sample& sample)
: AllScatterDevice(thisworld,sample) 
{
	p_sample->coordinate_sets.set_representation(CARTESIAN);
	
	string target = Params::Inst()->scattering.target;
	
	size_t NN = thisworld.size(); // Number of Nodes
	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();
    size_t NMYF = myframes.size();
	
    size_t NM = 1;
	if (Params::Inst()->scattering.average.orientation.vectors.size()>0) {
        NM = Params::Inst()->scattering.average.orientation.vectors.size();
    }
    size_t NMBLOCK = (NN<NM) ? NN : NM;
    
	if (thisworld.rank()==0) {	

        size_t memusage_scatmat = 2*sizeof(double)*NM*NMYF;

		Info::Inst()->write(string("Memory(Scattering Matrix): ")+to_s(memusage_scatmat)+string(" bytes"));

        // fault if not enough memory for scattering matrix
        if (memusage_scatmat>Params::Inst()->limits.memory.scattering_matrix) {
			Err::Inst()->write(string("Insufficient Buffer size for scattering matrix."));            
			Err::Inst()->write(string("Size required:")+to_s(memusage_scatmat)+string(" bytes"));            
			Err::Inst()->write(string("Configuration Parameter: limits.memory.scattering_matrix"));            
        }
	}	
}


void AllVectorsScatterDevice::scatter(size_t moffset,size_t mcount) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;
   
   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = myframes.size();
   size_t NM = mcount;
      
   p_a->resize(NM);
   for(size_t i = 0; i < NM; ++i)
   {
       (*p_a)[i].resize(NMYF,0);
   }
   
   string target = Params::Inst()->scattering.target;
   size_t NOA = p_sample->atoms.selections[target].indexes.size();
   
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       timer.start("sd:fs:f:ld");	
       CoordinateSet& cs = p_sample->coordinate_sets.load(myframes[fi]); 
       timer.stop("sd:fs:f:ld");	
           
  
       for(size_t mi = 0; mi < NM; ++mi)
       {
           CartesianCoor3D& q = qvectors[moffset+mi];

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
           (*p_a)[mi][fi] += complex<double>(Ar,Ai);
	   }
	}
}


void AllVectorsScatterDevice::init(CartesianCoor3D& q) {
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

size_t AllVectorsScatterDevice::get_numberofmoments() {
    return qvectors.size();
}

void AllVectorsScatterDevice::norm() {
    size_t NV = p_asingle->size();
    for(size_t i = 0; i < NV; ++i)
    {
        (*p_asingle)[i] /= qvectors.size();
    }
}

// end of file
