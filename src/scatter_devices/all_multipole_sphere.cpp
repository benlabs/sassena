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
#include "scatter_devices/all_multipole_sphere.hpp"

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


AllMSScatterDevice::AllMSScatterDevice(	boost::mpi::communicator scatter_comm,
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
	
	// init resolution 
	long LMAX = Params::Inst()->scattering.average.orientation.multipole.resolution;
    moments.push_back(make_pair<long,long>(0,0));
    for(long l = 1; l <= LMAX; ++l)
    {
        for(long m = -LMAX; m <= LMAX; ++m)
        {
            moments.push_back(make_pair<long,long>(l,m));
        }
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

void AllMSScatterDevice::scatter(size_t moffset,size_t mcount) {
   // outer loop: frames
   // inner loop: block of moments
   using namespace boost::numeric::ublas::detail;

   std::vector<double>& sfs = scatterfactors.get_all();

   size_t NMYF = myframes.size();
   size_t NM = mcount;
      
   double ql = q.length();
   double M_PI_four = 4*M_PI;
  
   string target = Params::Inst()->scattering.target;
   size_t NOA = p_sample->atoms.selections[target].indexes.size();
   
   p_a->assign(NM,std::vector<complex<double> >(NMYF,0));
   
   for(size_t fi = 0; fi < NMYF; ++fi)
   {
       timer.start("sd:fs:f:ld");	
       CoordinateSet& cs = *csets[fi]; 
       timer.stop("sd:fs:f:ld");	

       for(size_t mi = 0; mi < NM; ++mi)
       {          
           long l = moments[moffset+mi].first;
           long m = moments[moffset+mi].second;

           complex<double> A = 0;
	       for(size_t j = 0; j < NOA; ++j) {
	           double r   = cs.c1[j];
	 	       double phi = cs.c2[j];
		       double theta = cs.c3[j];
			
		       double esf = sfs[j];
               double p = ql*r;
           
			   complex<double> fmpiilesf = M_PI_four*pow(complex<double>(0,1.0),l) * esf;
			   double aabess = boost::math::sph_bessel(l,p);
		
			   complex<double> aa = conj(boost::math::spherical_harmonic(l,m,theta,phi)); 
               A += fmpiilesf * aabess* aa;
		   }
           (*p_a)[mi][fi] = A;

	   }
	}
}

size_t AllMSScatterDevice::get_numberofmoments() {
    return moments.size();	
}

void AllMSScatterDevice::norm() {
    size_t NV = p_asingle->size();
    for(size_t i = 0; i < NV; ++i)
    {
        (*p_asingle)[i] /= 4*M_PI;
    }
}

void AllMSScatterDevice::init(CartesianCoor3D& qin) {
    q = qin;
}

// end of file
