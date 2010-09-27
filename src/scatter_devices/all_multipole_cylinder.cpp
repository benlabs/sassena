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
#include "scatter_devices/all_multipole_cylinder.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>

// other headers
#include "math/coor3d.hpp"
#include "math/smath.hpp"
#include "decomposition/decompose.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;


AllMCScatterDevice::AllMCScatterDevice(	boost::mpi::communicator scatter_comm,
		boost::mpi::communicator fqt_comm,
		Sample& sample,
		vector<pair<size_t,CartesianCoor3D> > QIV,
		boost::asio::ip::tcp::endpoint fileserver_endpoint)
: AllScatterDevice(scatter_comm,fqt_comm,sample,QIV,fileserver_endpoint)
{
	p_sample->coordinate_sets.set_representation(SPHERICAL);
	
	string target = Params::Inst()->scattering.target;
	
	
	size_t NN = fqt_comm.size(); // Number of Nodes
	//size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	//size_t NF = sample.coordinate_sets.size();
    size_t NMYF = myframes.size();
	
	// init moments / resolution
    long LMAX = Params::Inst()->scattering.average.orientation.multipole.resolution;
    moments.push_back(0);
    for(long l = 1; l <= 4*LMAX; ++l)
    {
        moments.push_back(l);
    }
	
	size_t NM = moments.size();
    size_t NMBLOCK = (NN<NM) ? NN : NM;
    
	if (fqt_comm.rank()==0) {

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

void AllMCScatterDevice::scatter(size_t moffset,size_t mcount) {
   // outer loop: frames
   // inner loop: block of moments
   
   using namespace boost::numeric::ublas::detail;

   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = myframes.size();
   size_t NM = mcount;
   
   
   string target = Params::Inst()->scattering.target;
   size_t NOA = p_sample->atoms.selections[target].indexes.size();
   
   p_a->assign(NM,std::vector<complex<double> >(NMYF,0));

   CartesianCoor3D o = Params::Inst()->scattering.average.orientation.multipole.axis;
   o = o / o.length();

   // get the part of the scattering vector perpenticular to the o- orientation
   CartesianCoor3D qparallel = (o*q)*o; 
   CartesianCoor3D qperpenticular = q - qparallel; // define this as phi=0 
   double qr = qperpenticular.length();
   double qz = qparallel.length();
   //double ql = q.length();
   
   // frames = outer loop b/c of function call
   // moments = 2nd outer loop b/c of l lookup
   
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       timer.start("sd:fs:f:ld");	
       CoordinateSet& cs = *csets[fi]; 
       timer.stop("sd:fs:f:ld");	

	   for(size_t mi = 0; mi < NM; ++mi)
	   {
	       long l = moments[moffset+mi];
//           long l = (mi+3)/4; // back mapping mi -> l

           int func = 0; // mi+moffset == 0 -> A[0]
           if ((moffset+mi)!=0) {
               if ((mi%4)==1) func = 1; // A
               if ((mi%4)==2) func = 2; // B
               if ((mi%4)==3) func = 3; // C
               if ((mi%4)==0) func = 4; // D
           }

           // do an individual code block for each case of mi (A[0],A[x],B[x],C[x],D[x])
           // this avoids if statements within the fastest loop
           complex<double> A = 0;
           if (func==0) {
    	       for(size_t j = 0; j < NOA; ++j) {
    	           double r   = cs.c1[j];
    	 	       //double phi = cs.c2[j];
    		       double z   = cs.c3[j];

    		       double esf = sfs[j];

        	       //double psiphi = phi; // review this!

        	       double parallel_sign = 1.0;
        	       if ((z!=0) && (qz!=0)) {
        	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
        	       }

        	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
                   double p = r*qr;

                   // if func==0
            	   A += expi * (double)boost::math::cyl_bessel_j(0,p) * esf ;               
                }
           } else if (func==1) {
    	       for(size_t j = 0; j < NOA; ++j) {
    	           double r   = cs.c1[j];
    	 	       double phi = cs.c2[j];
    		       double z   = cs.c3[j];

    		       double esf = sfs[j];

        	       double psiphi = phi; // review this!

        	       double parallel_sign = 1.0;
        	       if ((z!=0) && (qz!=0)) {
        	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
        	       }

        	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
                   double p = r*qr;

                   // if func==1
            	   complex<double> fac1 = 2.0*powf(-1.0,l)*boost::math::cyl_bessel_j(2*l,p);
                   A += sqrt(0.5)*fac1 * expi *cos(2*l*psiphi) * esf;
       		                  
                }
           } else if (func==2) {
               for(size_t j = 0; j < NOA; ++j) {
    	           double r   = cs.c1[j];
    	 	       double phi = cs.c2[j];
    		       double z   = cs.c3[j];

    		       double esf = sfs[j];

        	       double psiphi = phi; // review this!

        	       double parallel_sign = 1.0;
        	       if ((z!=0) && (qz!=0)) {
        	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
        	       }

        	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
                   double p = r*qr;

                   // if func==2
            	   complex<double> fac1 = 2.0*powf(-1.0,l)*boost::math::cyl_bessel_j(2*l,p);
                   A += sqrt(0.5)*fac1 * expi *sin(2*l*psiphi) * esf;          
                }
           } else if (func==3) {
               for(size_t j = 0; j < NOA; ++j) {
    	           double r   = cs.c1[j];
    	 	       double phi = cs.c2[j];
    		       double z   = cs.c3[j];

    		       double esf = sfs[j];

        	       double psiphi = phi; // review this!

        	       double parallel_sign = 1.0;
        	       if ((z!=0) && (qz!=0)) {
        	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
        	       }

        	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
                   double p = r*qr;

                   // if func==3
            	   complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*boost::math::cyl_bessel_j(2*l-1,p));
                   A += sqrt(0.5)*fac2 * expi *cos((2*l-1)*psiphi) * esf; 
                }
           } else if (func==4) {
               for(size_t j = 0; j < NOA; ++j) {
    	           double r   = cs.c1[j];
    	 	       double phi = cs.c2[j];
    		       double z   = cs.c3[j];

    		       double esf = sfs[j];

        	       double psiphi = phi; // review this!

        	       double parallel_sign = 1.0;
        	       if ((z!=0) && (qz!=0)) {
        	           parallel_sign = (z*qz) / (abs(z)*abs(qz));			
        	       }

        	       complex<double> expi = exp(complex<double>(0,parallel_sign*z*qz));
                   double p = r*qr;

                   // if func==4
            	   complex<double> fac2 = complex<double>(0,1.0)*double(2.0*powf(-1.0,l-1)*boost::math::cyl_bessel_j(2*l-1,p));
                   A += sqrt(0.5)*fac2 * expi *sin((2*l-1)*psiphi) * esf;            
                }
            }
            (*p_a)[mi][fi] = A;
            
           } // mto
	} // NMYF
}

size_t AllMCScatterDevice::get_numberofmoments() {
    return moments.size();	
}

void AllMCScatterDevice::norm() {
    // do nothing
}

void AllMCScatterDevice::init(CartesianCoor3D& qin) {
    q = qin;
}

// end of file
