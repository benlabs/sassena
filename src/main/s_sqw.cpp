/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
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
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <fftw3.h>

// other headers
#include "math/coor3d.hpp"
#include "log.hpp"

#include "SassenaConfig.hpp"

using namespace std;

namespace po = boost::program_options;

void print_title() {
    
	Info::Inst()->write("This software is being developed by Benjamin Lindner.                    ");
	Info::Inst()->write("For help, suggestions or correspondense use:                             ");
	Info::Inst()->write("ben@benlabs.net, Benjamin Lindner (Main Developer, Impl. & Maintenance)  ");
	Info::Inst()->write("franc@cmm.ki.si, Franci Merzel (Methodology)                             ");
	Info::Inst()->write("For publications include the following references:                       ");
	Info::Inst()->write(".........................................................................");
	Info::Inst()->write("1. Sassena - Scattering Calculations on Parallel Computers               ");
	Info::Inst()->write("   to be published                                                       ");		
	Info::Inst()->write(".........................................................................");
    Info::Inst()->write(string("Version Information: ") + string(Sassena_VERSIONSTRING));
	Info::Inst()->write("");

}

void print_description() {
        Info::Inst()->write(".................................................................");
    	Info::Inst()->write("......................D.E.S.C.R.I.P.T.I.O.N......................");
    	Info::Inst()->write(".................................................................");	

    	Info::Inst()->write("This binary generates a trajectory index to allow random access");	
    	Info::Inst()->write("to frames within a molecular dynamics trajectory. ");
    	Info::Inst()->write(".................................................................");	
    	
}

bool init_commandline(int argc,char** argv,po::variables_map& vm) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("signal",po::value<string>()->default_value("signal.h5"),"name of the signal file")
        ("sqw",po::value<string>(),"name of the sqw file (defaults to signal file name with -sqw postfix)")
    ;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    
    if (vm.find("help")!=vm.end()) {
        cout << desc << endl;
        return false;
    }
    
    return true;
}

//void read_parameters(po::variables_map vm) {
//    Info::Inst()->write("Checking command line for parameter overwrite");
    // second stage , allow command line parameters to overwrite defaults in the configuration file
//}

void init_sqw(const string filename, size_t nf,size_t nq) {
    hid_t h5file = H5Fcreate( filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);   

    double fill_val_double = 0;

    size_t chunksize = 10000;

    {
        hsize_t dims[3];hsize_t mdims[3];hsize_t cdims[3];
        dims[0]=nq;              dims[1]=nf;     dims[2]=2;
        mdims[0]=H5S_UNLIMITED; mdims[1]=nf;     mdims[2]=2;
        size_t fqt_chunksize = chunksize;
        size_t fqt_chunksize2 = ((nf>0) && (nf<fqt_chunksize)) ? nf : fqt_chunksize;
        size_t fqt_chunksize1 = fqt_chunksize/fqt_chunksize2;
        cdims[0]=fqt_chunksize1; cdims[1]=fqt_chunksize2; cdims[2]=2;
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk( dcpl, 3, cdims);
        H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
        H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
        hid_t dspace = H5Screate_simple(3, dims, mdims); 
        hid_t ds = H5Dcreate(h5file, "sqw", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
        H5Pclose(lcpl);
        H5Pclose(dcpl);
        H5Pclose(dapl);
        H5Sclose(dspace);
        H5Dclose(ds);
    }   
    
    H5Fclose(h5file);
}

bool test_fqt(const std::string filename) {
        hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        
        if (h5file<0) {
            Err::Inst()->write("Error opening signal file. HDF5?");
            return false;
        }
        
        hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);

        if (ds_fqt<0) {
            Err::Inst()->write("Error opening signal file dataset 'fqt'");
            H5Fclose(h5file);            
            return false;
        }
        

        H5Dclose(ds_fqt);
        H5Fclose(h5file);

        return true;
}

size_t fqt_len(const std::string filename) {
        size_t retval =0;
        hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);    
        hid_t dspace_fqt = H5Dget_space(ds_fqt);

        hsize_t fqt_field[3];
        hsize_t fqt_field_max[3];

        H5Sget_simple_extent_dims(dspace_fqt,fqt_field,fqt_field_max);    

        retval= fqt_field[1];

        H5Sclose(dspace_fqt);    
        H5Dclose(ds_fqt);
        H5Fclose(h5file);

        return retval;
}

void write_frequencies(std::vector<double> frequencies, const std::string filename) {
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    
    hsize_t freq_field[1]; freq_field[0]=frequencies.size();
    H5LTmake_dataset_double(h5file,"frequencies",1,freq_field,&(frequencies[0]));

    H5Fclose(h5file);
}


std::vector<CartesianCoor3D> get_qvectors(std::string filename)
{
    herr_t status ;
    std::vector<CartesianCoor3D> qvectors;
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

    hsize_t dims[2];
    status = H5LTget_dataset_info(h5file,"qvectors",dims,NULL,NULL);
    std::vector<CartesianCoor3D> h5qvectors(dims[0]);
    if (dims[0]>0) {
        status = H5LTread_dataset_double(h5file,"qvectors",reinterpret_cast<double*>(&h5qvectors[0]));            
    }

    H5Fclose(h5file);

    return h5qvectors;
}

void write_qvectors(std::string filename,const std::vector<CartesianCoor3D>& qvectors) {
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    
    hsize_t qv_field[2]; 
    qv_field[0]=qvectors.size(); qv_field[1]=3;
    H5LTmake_dataset_double(h5file,"qvectors",2,qv_field,&(qvectors[0].x));

    H5Fclose(h5file);
}

void read_fqtslice(hid_t ds,size_t qindex,fftw_complex* workspace,size_t nf) {
    hid_t dspace = H5Dget_space(ds);

    hsize_t start[3];hsize_t count[3];
    hsize_t mdims[2];
    mdims[0]=nf; mdims[1]=2;
    hid_t mspace = H5Screate_simple(2, mdims, NULL);             
    
    start[0]=qindex;start[1]=0;start[2]=0;
    count[0]=1;count[1]=nf;count[2]=2;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
    start[0]=0;start[1]=0;
    count[0]=nf;count[1]=2;
    H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

    double* p_data = reinterpret_cast<double*>(workspace);
    H5Dread(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_data); 
        
    H5Sclose(dspace);    
}


void write_sqwslice(hid_t ds,size_t qindex, fftw_complex* workspace,size_t nf) {
    hid_t dspace = H5Dget_space(ds);
    
    hsize_t start[3];hsize_t count[3];

    hsize_t mdims[2];
    mdims[0]=nf; mdims[1]=2;
    hid_t mspace = H5Screate_simple(2, mdims, NULL);             

    
    start[0]=qindex;start[1]=0;start[2]=0;
    count[0]=1;count[1]=nf;count[2]=2;
    H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
    start[0]=0;start[1]=0;
    count[0]=nf;count[1]=2;
    H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

    double* p_data = reinterpret_cast<double*>(workspace);
    H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_data); 
        
    H5Sclose(dspace);    
}



//void write_sqw() {
//    size_t nq = data_fqt.size();
//    size_t nf = data_fqt[0]->size();
//    
//    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
//    
//    if (H5LTfind_dataset(h5file,"fqt")==1) {
//
//        double* p_fqtdata = (double*) malloc(sizeof(double)*2*nq*nf);
//        for(size_t i = 0; i < nq; ++i)
//        {
//            memcpy(&p_fqtdata[i*nf*2],&((*data_fqt[i])[0]),sizeof(double)*2*nf);
//        }
//        hsize_t dims[3];hsize_t start[3];hsize_t count[3];
//        H5LTget_dataset_info(h5file,"fqt",dims,NULL,NULL); 
//        
//        dims[0]+=nq;
//        hid_t ds = H5Dopen(h5file,"fqt",H5P_DEFAULT);
//        H5Dset_extent(ds,dims);
//        hid_t dspace = H5Dget_space(ds);
//
//        hsize_t mdims[3];
//        mdims[0]=nq; mdims[1]=nf;mdims[2]=2;
//        hid_t mspace = H5Screate_simple(3, mdims, NULL);             
//
//        start[0]=dims[0]-nq;start[1]=0;start[2]=0;
//        count[0]=nq;count[1]=nf;count[2]=2;
//        H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
//        start[0]=0;start[1]=0;start[2]=0;
//        count[0]=nq;count[1]=nf;count[2]=2;            
//        H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);
//
//        H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_fqtdata);  
//        
//        free(p_fqtdata);          
//    }
//}

int main(int argc,char* argv[]) {

	// make sure any singleton class exists:
	Info::Inst();
	Err::Inst();
	Warn::Inst();
	
	Info::Inst()->set_prefix(string("Info>>"));
	Warn::Inst()->set_prefix(string("Warn>>"));
	Err::Inst()->set_prefix(string("Err>>"));

    bool initstatus = true;
		
    po::variables_map vm;    
    print_title();
    print_description();
    initstatus = init_commandline(argc,argv,vm);
    if (!initstatus) {
        return 0;
    }

//    read_parameters(vm);

    using namespace boost::filesystem;
    path signalpath(vm["signal"].as<string>());
    if (!exists(signalpath)) {
        Err::Inst()->write(string("Signal file not found with path: ")+signalpath.filename().string());
        return -1;
    } else {
        Info::Inst()->write(string("Reading signal from file: ")+signalpath.filename().string());        
    }

    boost::filesystem::path sqwpath;
    if (vm.find("sqw")==vm.end()) {
        std::string tf = vm["signal"].as<string>();
        boost::filesystem::path fpath = tf;
    	string fdir;
        if (fpath.parent_path().is_complete()) {
    		fdir = fpath.parent_path().string();             
        } else {
    		fdir = (initial_path() / fpath).string();    	    
    	}
		sqwpath = (path(fdir).parent_path() / fpath.stem().string() ).string()+ string("-sqw.h5");
        Warn::Inst()->write("SQW file name not specified. Defaulting to signal filename with -sqw postfix:");
        Warn::Inst()->write(sqwpath.string());
    } else {
        sqwpath = vm["sqw"].as<string>();        
    }

    if (exists(sqwpath)) {
        int n=0;
        while (boost::filesystem::exists(sqwpath.filename().string()+".backup-"+boost::lexical_cast<string>(n))) n++;
        path newsqwpath = sqwpath.filename().string()+".backup-"+boost::lexical_cast<string>(n);
        Warn::Inst()->write(string("Moving old sqw file to ")+newsqwpath.filename().string());        
        boost::filesystem::rename(sqwpath,newsqwpath);
    }

    // try to open fqt in signal file
    Info::Inst()->write("Testing signal file");
    if (!test_fqt(signalpath.string())) {
        Err::Inst()->write("Error reading from signal file");
        throw;
    }

    size_t nq = 0;
    size_t nf = 0;
    size_t nfhalf = 0;
    // copy the qvectors to the sqw file
    {
        // gather size information, get nf
        Info::Inst()->write("Reading information from signal file");
        nf = fqt_len(signalpath.filename().string());
        nfhalf = nf/2;
        Info::Inst()->write(string("signal length nf=")+boost::lexical_cast<string>(nf));

        // reading qvectors, get nq
        Info::Inst()->write(string("Copying qvectors..."));    
        Info::Inst()->write(string("Reading from signal file"));    
        std::vector<CartesianCoor3D> qvectors = get_qvectors(signalpath.string());
        nq = qvectors.size();

        // init sqw, set nf, nq
        Info::Inst()->write("Initializing sqw file");    
        init_sqw(sqwpath.string(),nf,nq);
        
        // write qvectors
        Info::Inst()->write(string("Writing to sqw file"));    
        write_qvectors(sqwpath.string(),qvectors);
        
        // create frequency information from dt and nf
        Info::Inst()->write("Generating frequencies");
        double dt = 1;
        if (vm.find("dt")!=vm.end()) {
            dt = vm["dt"].as<double>();        
            if (dt==0) {
                Err::Inst()->write("dt must be larger than 0");                    
                throw;
            }
        } else {
            Warn::Inst()->write("No 'dt' specified. Assuming dt=1.");        
        }
        std::vector<double> frequencies(nf);
        for(size_t i = 0; i < nfhalf; ++i)
        {
            frequencies[i]=1.0*i/(2.0*nf*dt);
        }
        for(size_t i = nfhalf; i < nf; ++i)
        {
            frequencies[i]=-1.0*(nf-i)/(2.0*nf*dt);
        }

        // write frequencies
        Info::Inst()->write(string("Writing frequencies to sqw file"));    
        write_frequencies(frequencies,sqwpath.string());
    }

    std::vector<double> resolution(nf);
    {
        // default 10% gaussian
        double fwhm_t = 0.005*nf;
        double sigma_t = 2*fwhm_t / ( 2.0*sqrt(2.0*log(2.0)) );
        double mu = 0;
        for(size_t i = 0; i < nf; ++i)
        {
            resolution[i] = exp(-0.5*powf( (i - mu) / sigma_t,2));
        }        
    }

    {
            
        // compute for every q vector the fft for fqt
        hid_t fqtfile = H5Fopen( signalpath.string().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        hid_t sqwfile = H5Fopen( sqwpath.string().c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        
        hid_t ds_fqt = H5Dopen(fqtfile,"fqt",H5P_DEFAULT);
        hid_t ds_sqw = H5Dopen(sqwfile,"sqw",H5P_DEFAULT);
        
        // make workspace twice as large
        fftw_complex* workspace = (fftw_complex*) fftw_malloc((2*nf-1)*sizeof(fftw_complex)); 
        
        fftw_plan plan = fftw_plan_dft_1d((2*nf-1), workspace, workspace, FFTW_FORWARD, FFTW_ESTIMATE);    
        
        for(size_t i = 0; i < nq; ++i)
        {
            read_fqtslice(ds_fqt,i,workspace,nf);
            
            // apply resolution
            for(size_t n = 0; n < nf; ++n) {
                double a = workspace[n][0]*resolution[n];
                double b = workspace[n][1]*resolution[n];
                workspace[n][0]=a;
                workspace[n][1]=b;         
                workspace[2*nf-1-n][0]=a;      
                workspace[2*nf-1-n][1]=b;
            }
            
            fftw_execute_dft(plan,workspace,workspace);
                        
            write_sqwslice(ds_sqw,i,workspace,nf);
        }
        
        fftw_free(workspace);
        
        H5Dclose(ds_fqt);
        H5Dclose(ds_sqw);
        H5Fclose(fqtfile);
        H5Fclose(sqwfile);
    }


    
    // store the result into sqw file after each iteration

    Info::Inst()->write(string("Done."));
    
	return 0;
}


// end of file