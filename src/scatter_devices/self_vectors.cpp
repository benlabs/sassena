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
#include "scatter_devices/io/h5_fqt_interface.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"


using namespace std;


void SelfVectorsScatterDevice::compute() {

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

	m_spectrum.assign(NF,0);

    timer.start("sd:compute");
    // correlate or sum up
    if (Params::Inst()->scattering.correlation.type=="time") {

        for (long ai = 0; ai < m_indexes.size(); ai++) {
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
        	    for(size_t fi = 0; fi < NF; ++fi)
        	    {
        	    	m_spectrum[fi] += spectrum[fi];
        	    }
        	}
        }
    } else if (Params::Inst()->scattering.correlation.type=="infinite-time") {

        for (long ai = 0; ai < m_indexes.size(); ai++) {
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
                infinite_correlate();
    	        timer.stop("sd:correlate");

                vector<complex<double> >& spectrum = (*p_asingle);
        	    for(size_t fi = 0; fi < NF; ++fi)
        	    {
        	    	m_spectrum[fi] += spectrum[fi];
        	    }
        	}
    	}
    } else if (Params::Inst()->scattering.correlation.type=="average-time") {

        for (long ai = 0; ai < m_indexes.size(); ai++) {
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
                average_correlate();
    	        timer.stop("sd:correlate");

                vector<complex<double> >& spectrum = (*p_asingle);
        	    for(size_t fi = 0; fi < NF; ++fi)
        	    {
        	    	m_spectrum[fi] += spectrum[fi];
        	    }
        	}
    	}
    } else if (Params::Inst()->scattering.correlation.type=="instant-time") { // this corresponds to "static", i.e. instant correlation = conjmultiply, only elastic part of the function
                    // if not time correlated, the conjmultiply negates phase information
                    // this simplifies formulas

            //        p_asingle->assign(NF,0);

                    for (long ai = 0; ai < m_indexes.size(); ai++) {
                    	double s = scatterfactors.get(m_indexes[ai]);
                        for(size_t fi = 0; fi < NF; ++fi)
                        {
                            m_spectrum[fi] += NM*s*s;
                        }
                    }
    } else {
		Err::Inst()->write("Correlation type not understood.");        
        throw;
    }
    timer.stop("sd:compute");


	// these functions operate on m_spectrum:

    timer.start("sd:norm");
    norm();
    timer.stop("sd:norm");

	// finally assemble results on the head node:
    timer.start("sd:gather_sum");
    gather_sum(); // head node has everything in a (vector< complex<double> >)
    timer.stop("sd:gather_sum");

    m_writeflag = true;
}

void SelfVectorsScatterDevice::write() {
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

void SelfVectorsScatterDevice::next() {
	if (m_current_qvector>=m_qvectorindexpairs.size()) return;

	m_current_qvector+=1;
}



double SelfVectorsScatterDevice::progress() {
	size_t total = m_qvectorindexpairs.size();
	size_t current = m_current_qvector;
	double progress = 0.0;
	if (total==current) {
		progress = 1.0;
	} else if (current!=0) {
		progress = ((current*1.0)/total);
	}

	return progress;
}

SelfVectorsScatterDevice::SelfVectorsScatterDevice(
		boost::mpi::communicator scatter_comm,
		boost::mpi::communicator fqt_comm,
		Sample& sample,
		std::vector<std::pair<size_t,CartesianCoor3D> > QIV,
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

	p_sample->coordinate_sets.set_representation(CARTESIAN);	
	p_sample->coordinate_sets.set_selection(sample.atoms.selections[target]);


    // if no qvectors are to be calculated , DO not waste any time on transposing and/or loading the data..
    if (QIV.size()==0) {
        return;
    }

	p_asingle = new vector< complex<double> >; // initialize with no siz
	
	scatterfactors.set_sample(sample);
	scatterfactors.set_selection(sample.atoms.selections[target]);
	scatterfactors.set_background(true);


	size_t rank = m_fqtcomm.rank();
	size_t NN = m_fqtcomm.size(); // Number of Nodes

	size_t NA = sample.atoms.selections[target].indexes.size(); // Number of Atoms
	size_t NF = sample.coordinate_sets.size();

	EvenDecompose framesdecomp(NF,NN);
    vector<size_t> myframes = framesdecomp.indexes_for(rank);
    size_t NMYF = myframes.size();

	EvenDecompose indexdecomp(NA,NN);
    m_indexes = indexdecomp.indexes_for(rank);
    size_t NMYI = m_indexes.size();
    
	if (NN>NA) {
		Err::Inst()->write("Not enough atoms for decomposition");
		throw;
	}


    m_x.assign(NMYI,vector<double>(NF));
    m_y.assign(NMYI,vector<double>(NF));
    m_z.assign(NMYI,vector<double>(NF));


    
    // first do a global communication to distribute the atoms
    size_t seglength = indexdecomp.indexes_for(rank).size();
    size_t segoffset = indexdecomp.indexes_for(rank)[0];

    vector<size_t> slengths(NN);
    vector<size_t> soffsets(NN);

    boost::mpi::all_gather(m_fqtcomm,seglength,&(slengths[0]));
    boost::mpi::all_gather(m_fqtcomm,segoffset,&(soffsets[0]));
    // now slengths and soffsets index atom stripes

    // load all frames first

    timer.start("sd:data");
    vector<vector<double> > localdatax(NN);
    vector<vector<double> > localdatay(NN);
    vector<vector<double> > localdataz(NN);
    for(size_t i = 0; i < NN; ++i)
    {
        localdatax[i].resize(slengths[i]*NMYF);
        localdatay[i].resize(slengths[i]*NMYF);
        localdataz[i].resize(slengths[i]*NMYF);
    }

    for(size_t fi = 0; fi < NMYF; ++fi)
    {
        CoordinateSet* p_cset = p_sample->coordinate_sets.load(myframes[fi]);
	    for(size_t i = 0; i < NN; ++i)
	    {
            memcpy(&localdatax[i][fi*slengths[i]],&(p_cset->c1[soffsets[i]]),sizeof(double)*slengths[i]);
            memcpy(&localdatay[i][fi*slengths[i]],&(p_cset->c2[soffsets[i]]),sizeof(double)*slengths[i]);
            memcpy(&localdataz[i][fi*slengths[i]],&(p_cset->c3[soffsets[i]]),sizeof(double)*slengths[i]);	    	       
	    }
        delete p_cset;
    }

    timer.stop("sd:data");    

    timer.start("sd:data:trans");
    
    // first do a global communication to distribute the atoms
    size_t flength = framesdecomp.indexes_for(rank).size();
    size_t foffset = framesdecomp.indexes_for(rank)[0];

    vector<size_t> flengths(NN);
    vector<size_t> foffsets(NN);

    boost::mpi::all_gather(m_fqtcomm,flength,&(flengths[0]));
    boost::mpi::all_gather(m_fqtcomm,foffset,&(foffsets[0]));
    // now slengths and soffsets index atom stripes

    // next do a global communication to distribute the datasizes
    
    // now exchange data.
    for(size_t nodeiter = 0; nodeiter < NN; ++nodeiter)
    {
        size_t thisnode = rank;
        size_t sendnode = (thisnode+nodeiter)%NN;
        size_t recvnode = ((thisnode+NN)-nodeiter)%NN;

        if (nodeiter==0) {
            
            for(size_t i = 0; i < flengths[thisnode]; ++i)
            {
                size_t abs_framepos = foffsets[recvnode]+i;
                size_t rel_framepos_offset = i*slengths[thisnode];
                for(size_t j = 0; j < slengths[thisnode]; ++j)
                {
                    size_t rel_framepos = rel_framepos_offset + j;
                    m_x[j][abs_framepos]=localdatax[thisnode][rel_framepos];
                    m_y[j][abs_framepos]=localdatay[thisnode][rel_framepos];
                    m_z[j][abs_framepos]=localdataz[thisnode][rel_framepos];
                }
            }
            
        } else {
            std::list< boost::mpi::request > requests;
            
            vector<double> recvxcoords(flengths[recvnode]*slengths[thisnode]);
            vector<double> recvycoords(flengths[recvnode]*slengths[thisnode]);
            vector<double> recvzcoords(flengths[recvnode]*slengths[thisnode]);

            // non-blocking recv
            requests.push_back(m_fqtcomm.irecv(recvnode,0,&(recvxcoords[0]),flengths[recvnode]*slengths[thisnode]));
            requests.push_back(m_fqtcomm.irecv(recvnode,0,&(recvycoords[0]),flengths[recvnode]*slengths[thisnode]));
            requests.push_back(m_fqtcomm.irecv(recvnode,0,&(recvzcoords[0]),flengths[recvnode]*slengths[thisnode]));
           
            // blocking send
            double* p_xsegment = (double*) &(localdatax[sendnode][0]);
            double* p_ysegment = (double*) &(localdatay[sendnode][0]);
            double* p_zsegment = (double*) &(localdataz[sendnode][0]);
            
            m_fqtcomm.send(sendnode,0,p_xsegment,flengths[thisnode]*slengths[sendnode]);
            m_fqtcomm.send(sendnode,0,p_ysegment,flengths[thisnode]*slengths[sendnode]);
            m_fqtcomm.send(sendnode,0,p_zsegment,flengths[thisnode]*slengths[sendnode]);
            
            boost::mpi::wait_all(requests.begin(),requests.end());

            for(size_t i = 0; i < flengths[recvnode]; ++i)
            {
                size_t abs_framepos = foffsets[recvnode]+i;
                size_t rel_framepos_offset = i*slengths[thisnode];
                for(size_t j = 0; j < slengths[thisnode]; ++j)
                {
                    size_t rel_framepos = rel_framepos_offset + j;
                    m_x[j][abs_framepos]=recvxcoords[rel_framepos];
                    m_y[j][abs_framepos]=recvycoords[rel_framepos];
                    m_z[j][abs_framepos]=recvzcoords[rel_framepos];
                }
            }
        }

    }
    timer.stop("sd:data:trans");
}

void SelfVectorsScatterDevice::scatter(size_t ai, size_t mi) {

    vector<double>& x = m_x[ai];
    vector<double>& y = m_y[ai];
    vector<double>& z = m_z[ai];
    
	// this is broken <-- revise this!!!
	double s = scatterfactors.get(m_indexes[ai]);
    
    size_t NF = p_sample->coordinate_sets.size();
	
    p_asingle->assign(NF,0);	

    double qx = qvectors[mi].x;
    double qy = qvectors[mi].y;
    double qz = qvectors[mi].z;
    
	for(size_t j = 0; j < NF; ++j)
	{
		double x1 = x[j];
		double y1 = y[j];
		double z1 = z[j];	

		double p1 = x1*qx+y1*qy+z1*qz;
        
        double sp1 = sin(p1);
        double cp1 = cos(p1);
        (*p_asingle)[j] = s*complex<double>(cp1,sp1);
	}
	
}

void SelfVectorsScatterDevice::multiply_alignmentfactors(size_t mi) {
    
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
    
    std::vector< std::complex<double> >* p_correlated_a = new std::vector< std::complex<double> >(NF);
      
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


void SelfVectorsScatterDevice::average_correlate() {
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

void SelfVectorsScatterDevice::infinite_correlate() {
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

void SelfVectorsScatterDevice::gather_sum() {
    size_t NN = m_fqtcomm.size();
    size_t NF = p_sample->coordinate_sets.size();
        
    if (NF<1) return;
    
    if (m_fqtcomm.rank()==0) {
        vector< complex<double> > received_a(NF);
        double* p_a_double = (double*) &(received_a.at(0));   
        vector< complex<double> >& a = m_spectrum;
        
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
        double* p_a_double = (double*) &(m_spectrum[0]);   
        m_fqtcomm.send(0,0,p_a_double,2*NF);
        m_spectrum.clear();
    }

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
