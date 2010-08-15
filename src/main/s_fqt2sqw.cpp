/*
 *  sassena.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "common.hpp"

// standard header
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/regex.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/math/distributions/normal.hpp>
#include <fftw3.h>
#include <H5Cpp.h>

// other headers
#include "log.hpp"
#include "scatter_devices/scatter_spectrum.hpp"

using namespace std;

void init_output(std::string filename, vector<size_t> qindexes,vector<CartesianCoor3D> qvectors) {
    
    H5::H5File sqwfile(filename,H5F_ACC_TRUNC);
    hsize_t dims1[2],maxdims1[2],cdims1[2];
    dims1[0]=qindexes.size();
    dims1[1]=3;
    maxdims1[0]=H5S_UNLIMITED;
    maxdims1[1]=3;
    H5::DSetCreatPropList qvector_cparms;
    cdims1[0]=1;
    cdims1[1]=1;
    qvector_cparms.setChunk(2,cdims1);
    qvector_cparms.setFillValue( H5::PredType::NATIVE_DOUBLE, &fill_val_double);
    sqwfile.createDataSet( "qvectors", H5::DataType(H5::PredType::NATIVE_DOUBLE), H5::DataSpace(2,dims1,maxdims1),qvector_cparms );

    hsize_t dims2[3],maxdims2[3],cdims2[3];
    dims2[0]=qvectors.size();
    dims2[1]=nf;
    dims2[2]=2;
    maxdims2[0]=H5S_UNLIMITED;
    maxdims2[1]=nf;
    maxdims2[2]=2;
    H5::DSetCreatPropList sqw_cparms;
    cdims2[0]=1;
    cdims2[1]=1;
    cdims2[2]=1;
    sqw_cparms.setChunk(3,cdims2);
    sqw_cparms.setFillValue( H5::PredType::NATIVE_DOUBLE, &fill_val_double);
    sqwfile.createDataSet( "sqw", H5::DataType(H5::PredType::NATIVE_DOUBLE), H5::DataSpace(3,dims2,maxdims2),sqw_cparms );

    hsize_t foffset1[2];  // Start of hyperslab
    foffset1[0]=0;
    foffset1[1]=0;
    fspace = ds_qv.getSpace();
    fspace.selectHyperslab(H5S_SELECT_SET,dims1,foffset1);
    mspace = H5::DataSpace(2,dims1);
    ds_qv.write(reinterpret_cast<double*>(const_cast<CartesianCoor3D*>(&qvectors[0])), H5::PredType::NATIVE_DOUBLE,mspace,fspace);
    sqwfile.close();
}

void store_output(size_t qindex,vector<complex<double> >& sqw) {
    
}

int main(int argc,char** argv) {

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(string("Info>>"));
	Warn::Inst()->set_prefix(string("Warn>>"));
	 Err::Inst()->set_prefix(string("Err>>"));
	
     po::options_description desc("Allowed options");
     desc.add_options()
         ("help", "produce help message")
         ("input", po::value<string>()->default_value("fqt.h5"),  "name of the data input file")
         ("output",po::value<string>()->default_value("sqw.h5"),"name of the data output file")
     ;

     po::variables_map vm;
     po::store(po::parse_command_line(argc, argv, desc), vm);
     po::notify(vm);

     if (vm.count("help")) {
         cout << desc << endl;
         return 1;
     }

     if (vm.count("input")==0) {
         Info::Inst()->write("Require input file");
         cout << desc << endl;
     }
     
     H5::H5File h5file;

     h5file.openFile(vm["input"].as<string>(),H5F_ACC_RDONLY);

     H5::DataSet ds_qv = h5file.openDataSet("qvectors");
     H5::DataSet ds_okstatus = h5file.openDataSet("okstatus");
     H5::DataSet ds_fqt = h5file.openDataSet("fqt");

     hsize_t qvector_field[2];
     ds_qv.getSpace().getSimpleExtentDims(qvector_field);

     std::vector<CartesianCoor3D> h5qvectors(qvector_field[0]);
     ds_qv.read(reinterpret_cast<double*>(&h5qvectors[0]), H5::PredType::NATIVE_DOUBLE);

     hsize_t okstatus_field[1];
     ds_okstatus.getSpace().getSimpleExtentDims(okstatus_field);
     
     std::vector<int> okstatus(okstatus_field[0]);
     ds_okstatus.read(reinterpret_cast<int*>(&okstatus[0]), H5::PredType::NATIVE_INT);
     
     hsize_t fqt_field[3];
     ds_fqt.getSpace().getSimpleExtentDims(fqt_field);
     
     
     vector<size_t> qindexes;
     for(size_t i = 0; i < okstatus.size(); ++i)
     {
         if (okstatus[i]==1) {
             qindexes.push_back(i);
         }
     }

     init_output(vm["output"].as<string>(),qindexes,h5qvectors);

     for(size_t i = 0; i < qindexes.size(); ++i)
     {
        
     }
     
    
    ss.read_plain(argv[1],"txt");
    ss.write_plain("test.h5","hdf5");
    return 0;
    if (ss.size()<1) return 0;
    
    size_t seqlen = ss.size();
    
    // check sequence length
    for(size_t si = 0; si < ss.size(); ++si)
    {
        if (seqlen!=ss.amplitudes[si].size()) {
            Err::Inst()->write("Different Time Sequence Lengths. Exiting");
            return 1;
        }
    }
    
    size_t NF = seqlen;
    
 //   double dE = boost::lexical_cast<double>(argv[3]);
     double sigma_t = boost::lexical_cast<double>(argv[3]);
 //    double sigma_e = 4*dE/(3*M_PI);
 //    double sigma_t = 2/(M_PI*sigma_e);
//    double sigma_e = dE/2.354820045;
//    double sigma_t = 1.0/(1.5192669*sigma_e);
    boost::math::normal normdist(0,sigma_t);
    vector<double> rt;

    cout << "Generating Gaussian" << endl;
    // expand rt on demand
    if (NF>rt.size()) {
        size_t rtsize = rt.size();
        for(size_t i = rtsize; i < NF-rtsize; ++i) { 
            rt.push_back(boost::math::pdf(normdist,i)*sqrt(2*M_PI)*sigma_t);
        }
    }


    bool periodic = boost::lexical_cast<bool>(argv[4]);

    cout << "Fourier Transform" << endl;

    if (periodic) {

        fftw_complex *wspace;
        fftw_plan p1;
        wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NF-2));
        p1 = fftw_plan_dft_1d(2*NF-2, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
    
        for(size_t si = 0; si < ss.size(); ++si)
        {
            cout << "QVector " << si << endl;

            Info::Inst()->write(string("qvector ")+to_s(si));
            vector<complex<double> >& fqt = ss.amplitudes[si];
                                    
            for(size_t i = 0; i < NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[i].real()*rt[i]; 
                wspace[i][1]= fqt[i].imag()*rt[i];                            
            }

            for(size_t i = NF; i < 2*NF-2; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[2*NF-2-i].real()*rt[2*NF-2-i]; 
                wspace[i][1]= fqt[2*NF-2-i].imag()*rt[2*NF-2-i];                            
            }

            fftw_execute(p1); /* repeat as needed */
            for(size_t i = 0; i < NF; ++i) {
                fqt[i]=0.5*complex<double>(wspace[i][0],wspace[i][1])*( 1.0 / ( 2*NF-2 ) );
            }
        }
        
        fftw_destroy_plan(p1);
        fftw_free(wspace);

    } else {
        fftw_complex *wspace;
        fftw_plan p1;
        wspace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NF));
        p1 = fftw_plan_dft_1d(2*NF, wspace, wspace, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for(size_t si = 0; si < ss.size(); ++si)
        {
            Info::Inst()->write(string("qvector ")+to_s(si));
            vector<complex<double> >& fqt = ss.amplitudes[si];
            
            for(size_t i = 0; i < NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= fqt[i].real()*rt[i]; 
                wspace[i][1]= fqt[i].imag()*rt[i];                            

            }

            for(size_t i = NF; i < 2*NF; ++i) {
                // multiply w/ gaussian here = resolution function
                wspace[i][0]= 0; 
                wspace[i][1]= 0;                            

            }
            fftw_execute(p1); /* repeat as needed */
            for(size_t i = 0; i < NF; ++i) {
                fqt[i]=0.5*complex<double>(wspace[i][0],wspace[i][1])*( 1.0 / ( 2*NF ) );
            }

        }        
        
        fftw_destroy_plan(p1);
        fftw_free(wspace);

    }

    cout << "Writing SQW: " << argv[2] << endl;

    std::ofstream ofile(argv[2]);
    
    for(size_t si = 0; si < ss.size(); ++si) {
        vector<complex<double> >& fqt = ss.amplitudes[si];
        for(size_t i = 0; i < fqt.size(); ++i)
        {
            ofile << fqt[i].real() << endl;
        }
    }
    cout << "done" << endl;
    

	return 0;
}

// end of file