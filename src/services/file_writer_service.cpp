/*
 *  file_writer_service.cpp
 *
 *  Created on: Sep 04, 2009
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "services/file_writer_service.hpp"
#include <assert.h>


#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <log.hpp>
#include <control.hpp>
#include <math/smath.hpp>
using namespace std;

std::vector<size_t> HDF5WriterService::init(const std::vector<CartesianCoor3D>& qvectors,size_t nf) {

if (nf<1) {
    Err::Inst()->write("Number of frames must not be negative or zero!");    
    throw;
}


if (boost::filesystem::exists(m_filename)) {
    htri_t h5trialcode = H5Fis_hdf5(m_filename.c_str());  
    if (h5trialcode>0) {
        return init_reuse(qvectors,nf);
     } else if (h5trialcode==0) {
         Err::Inst()->write(m_filename + string(" does not to be a HDF5 data file. aborting..."));
         throw;
     } else {
         Err::Inst()->write("Error when testing for HDF5 character of data file. aborting...");
         throw;         
     }
    } else {
        return init_new(qvectors,nf);
    }
}



std::vector<size_t> HDF5WriterService::init_new(const std::vector<CartesianCoor3D>& qvectors,size_t nf)
{     
    const string filename = m_filename;
    hid_t h5file = H5Fcreate( filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);   

    int fill_val_int = 0;
    double fill_val_double = 0;
    
    // qvector
    
    hsize_t dims1[2],maxdims1[2],cdims1[2];
    dims1[0]=qvectors.size();
    dims1[1]=3;
    maxdims1[0]=H5S_UNLIMITED;
    maxdims1[1]=3;
    cdims1[0]=10000;
    cdims1[1]=3;
    hid_t qvector_lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t qvector_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t qvector_dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( qvector_dcpl, 2, cdims1);
    H5Pset_fill_value( qvector_dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
    if (Params::Inst()->limits.data.alloc_early) {
        H5Pset_alloc_time(qvector_dcpl,H5D_ALLOC_TIME_EARLY);
    } else {
        H5Pset_alloc_time(qvector_dcpl,H5D_ALLOC_TIME_DEFAULT);
    }
    hid_t dspace_qv = H5Screate_simple(2, dims1, maxdims1); 
    hid_t ds_qv = H5Dcreate(h5file, "qvectors", H5T_NATIVE_DOUBLE, dspace_qv, qvector_lcpl,qvector_dcpl,qvector_dapl);
    H5Pclose(qvector_lcpl);
    H5Pclose(qvector_dcpl);
    H5Pclose(qvector_dapl);
    H5Sclose(dspace_qv);
    H5Dclose(ds_qv);

    // fqt

    hsize_t dims2[3],maxdims2[3],cdims2[3];
    dims2[0]=qvectors.size();
    dims2[1]=nf;
    dims2[2]=2;
    maxdims2[0]=H5S_UNLIMITED;
    maxdims2[1]=nf;
    maxdims2[2]=2;
    cdims2[0]=1;
    cdims2[1]= ((nf>0) && (nf<10000)) ? nf : 10000;
    cdims2[2]=2;
    hid_t fqt_lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t fqt_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t fqt_dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( fqt_dcpl, 3, cdims2);
    H5Pset_fill_value( fqt_dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
    if (Params::Inst()->limits.data.alloc_early) {
        H5Pset_alloc_time(fqt_dcpl,H5D_ALLOC_TIME_EARLY);
    } else {
        H5Pset_alloc_time(fqt_dcpl,H5D_ALLOC_TIME_DEFAULT);
    }
        

    hid_t dspace_fqt = H5Screate_simple(3, dims2, maxdims2); 
    hid_t ds_fqt = H5Dcreate(h5file, "fqt", H5T_NATIVE_DOUBLE, dspace_fqt, fqt_lcpl,fqt_dcpl,fqt_dapl);
    H5Pclose(fqt_lcpl);
    H5Pclose(fqt_dcpl);
    H5Pclose(fqt_dapl);
    H5Sclose(dspace_fqt);
    H5Dclose(ds_fqt);
    
    // fq
    
    
    hsize_t dims_fq[2],maxdims_fq[2],cdims_fq[2];
    dims_fq[0]=qvectors.size();
    dims_fq[1]=2;
    maxdims_fq[0]=H5S_UNLIMITED;
    maxdims_fq[1]=2;
    cdims_fq[0]= 10000;
    cdims_fq[1]=2;
    hid_t fq_lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t fq_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t fq_dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( fq_dcpl, 2, cdims_fq);
    H5Pset_fill_value( fq_dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
    if (Params::Inst()->limits.data.alloc_early) {
        H5Pset_alloc_time(fq_dcpl,H5D_ALLOC_TIME_EARLY);
    } else {
        H5Pset_alloc_time(fq_dcpl,H5D_ALLOC_TIME_DEFAULT);
    }
    hid_t dspace_fq = H5Screate_simple(2, dims_fq, maxdims_fq); 
    hid_t ds_fq = H5Dcreate(h5file, "fq", H5T_NATIVE_DOUBLE, dspace_fq, fq_lcpl,fq_dcpl,fq_dapl);
    H5Pclose(fq_lcpl);
    H5Pclose(fq_dcpl);
    H5Pclose(fq_dapl);
    H5Sclose(dspace_fq);
    H5Dclose(ds_fq);

    
    // checkpoint
    
    hsize_t dims3[1],maxdims3[1],cdims3[1];
    dims3[0]=qvectors.size();
    maxdims3[0]=H5S_UNLIMITED;
    cdims3[0]=10000;
    hid_t cp_lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t cp_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t cp_dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( cp_dcpl, 1, cdims3);
    H5Pset_fill_value( cp_dcpl, H5T_NATIVE_INT, &fill_val_int);
    if (Params::Inst()->limits.data.alloc_early) {
        H5Pset_alloc_time(cp_dcpl,H5D_ALLOC_TIME_EARLY);
    } else {
        H5Pset_alloc_time(cp_dcpl,H5D_ALLOC_TIME_DEFAULT);
    }
    hid_t dspace_checkpoint = H5Screate_simple(1, dims3, maxdims3); 
    hid_t ds_checkpoint = H5Dcreate(h5file, "checkpoint", H5T_NATIVE_INT, dspace_checkpoint, cp_lcpl,cp_dcpl,cp_dapl);
    H5Pclose(cp_lcpl);
    H5Pclose(cp_dcpl);
    H5Pclose(cp_dapl);
    H5Sclose(dspace_checkpoint);
    H5Dclose(ds_checkpoint);
 
    std::vector<size_t> qindexes;   
    for(size_t i = 0; i < qvectors.size(); ++i) qindexes.push_back(i);

    // (over)write qvectors values for index positions
    hsize_t foffset1[2];
    foffset1[0]=0;
    foffset1[1]=0;
    hid_t ds_qv2 = H5Dopen(h5file,"qvectors",H5P_DEFAULT);
    hid_t fspace = H5Dget_space(ds_qv2);
    H5Sselect_hyperslab(fspace,H5S_SELECT_SET,foffset1,NULL,dims1,NULL);
    H5Dwrite(ds_qv2,H5T_NATIVE_DOUBLE,H5S_ALL,fspace,H5P_DEFAULT,reinterpret_cast<double*>(const_cast<CartesianCoor3D*>(&qvectors[0])));
    
    H5Sclose(fspace);
    H5Dclose(ds_qv2);    
    H5Fclose(h5file);
    
    return qindexes;
}

bool HDF5WriterService::test_fqt_dim( size_t nf) {
    const std::string filename = m_filename;
    bool retval = false;
    //////////////////////////////
    // read info from data file
    //////////////////////////////
    hid_t h5file = H5Fopen( m_filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);
    hid_t dspace_fqt = H5Dget_space(ds_fqt);
    
    hsize_t fqt_field[3];
    hsize_t fqt_field_max[3];
    
    H5Sget_simple_extent_dims(dspace_fqt,fqt_field,fqt_field_max);    
    if (fqt_field[1]==nf) retval=true;    

    H5Sclose(dspace_fqt);    
    H5Dclose(ds_fqt);
    H5Fclose(h5file);
    
    return retval;
}


std::vector<size_t> HDF5WriterService::init_reuse(const std::vector<CartesianCoor3D>& qvectors,size_t nf)
{     
    const std::string filename = m_filename;
    Warn::Inst()->write("Data file re-use currently broken. will create a new file.");
    if (true) {
//    if (test_fqt_dim(nf)==false) {
        // this in an "incompatible" data file
        // mv old one to a backup
        int n=0;
        while (boost::filesystem::exists(filename+".backup-"+boost::lexical_cast<string>(n))) n++;
        std::string newfilename = filename+".backup-"+boost::lexical_cast<string>(n);
        Warn::Inst()->write(string("Moving old data file to ")+newfilename);        
        boost::filesystem::rename(filename,newfilename);
        // escape by returning init_new
        return init_new(qvectors,nf);  
    }
    
    //////////////////////////////
    // read info from data file
    //////////////////////////////
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t ds_qv = H5Dopen(h5file,"qvectors",H5P_DEFAULT);
    hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);
    hid_t ds_fq = H5Dopen(h5file,"fq",H5P_DEFAULT);
    hid_t ds_checkpoint = H5Dopen(h5file,"checkpoint",H5P_DEFAULT);
        
    hid_t dspace_qv = H5Dget_space(ds_qv);
    hid_t dspace_fqt = H5Dget_space(ds_fqt);
    hid_t dspace_fq = H5Dget_space(ds_fq);
    hid_t dspace_checkpoint = H5Dget_space(ds_checkpoint);
    
    hsize_t qvector_field[2];
    hsize_t checkpoint_field[1];
    hsize_t fqt_field[3];
    hsize_t fq_field[2];
    hsize_t qvector_field_max[2];
    hsize_t checkpoint_field_max[1];
    hsize_t fqt_field_max[3];
    hsize_t fq_field_max[2];
    
    H5Sget_simple_extent_dims(dspace_qv,qvector_field,qvector_field_max);
    H5Sget_simple_extent_dims(dspace_fqt,fqt_field,fqt_field_max);
    H5Sget_simple_extent_dims(dspace_fq,fq_field,fq_field_max);
    H5Sget_simple_extent_dims(dspace_checkpoint,checkpoint_field,checkpoint_field_max);
    
    H5Sclose(dspace_qv);
    H5Sclose(dspace_fqt);
    H5Sclose(dspace_fq);
    H5Sclose(dspace_checkpoint);
    
    std::vector<CartesianCoor3D> h5qvectors(qvector_field[0]);
    H5Dread(ds_qv,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,reinterpret_cast<double*>(&h5qvectors[0]));
    
    std::vector<int> checkpoint(checkpoint_field[0]);
    H5Dread(ds_checkpoint,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,reinterpret_cast<int*>(&checkpoint[0]));

    //////////////////////////////
    // detect faulty entries
    // re-use those entries which have a 0 at the checkpoint
    //////////////////////////////
    std::vector<size_t> qindexes;
    std::vector<CartesianCoor3D> oldqvectors;
    
    for(size_t i = 0; i < checkpoint.size(); ++i)
    {
        if (checkpoint[i]==0) {
            qindexes.push_back(i);                 
        } else {
            oldqvectors.push_back(h5qvectors[i]);
        }
    }

    std::vector<CartesianCoor3D> finalqvectors;

    //////////////////////////////
    // build a list of new qvectors.
    // disregard those which are already computed (oldqvectors)
    //////////////////////////////
    
    for(size_t i = 0; i < qvectors.size(); ++i)
    {
        bool oldqfound=false;
        for(size_t j = 0; j< oldqvectors.size(); ++j)
        {
            if (fabs(qvectors[i].x-oldqvectors[j].x)>10*numeric_limits<double>::min()) continue;
            if (fabs(qvectors[i].y-oldqvectors[j].y)>10*numeric_limits<double>::min()) continue;
            if (fabs(qvectors[i].z-oldqvectors[j].z)>10*numeric_limits<double>::min()) continue;
            oldqfound=true; break;
        }
        if (!oldqfound) finalqvectors.push_back(qvectors[i]);
    }
    
    //////////////////////////////
    // finalqvectors contains only new qvectors
    // qindexes contains list of possible storage points
    // 3 possiblities:
    // size(finalqvectors)>size(qindexes). need to extend data sets
    // size(finalqvectors)==size(qindexes). No extension necessary
    // size(finalqvectors)<size(qindexes). reduce qindexes to size of finalqvectors
    //////////////////////////////

    if (finalqvectors.size()>0) {
        size_t extendsize=0;
        if (finalqvectors.size()<qindexes.size()) {
            qindexes.resize(finalqvectors.size());
        } else {
            extendsize = finalqvectors.size()-qindexes.size();                
        }

        if (extendsize>0) {
            for(size_t i = 0; i < extendsize; ++i)
            {
                qindexes.push_back(qvector_field[0]+i);
            }
            qvector_field[0]+=extendsize;
            H5Dset_extent(ds_qv,qvector_field);

            checkpoint_field[0]+=extendsize;
            H5Dset_extent(ds_checkpoint,checkpoint_field);

            fqt_field[0]+=extendsize;
            H5Dset_extent(ds_fqt,fqt_field);
            
            fq_field[0]+=extendsize;
            H5Dset_extent(ds_fq,fq_field);
        }    

        //////////////////////////////
        // qindexes contain the indexes for the qvector storage points
        // finalqvectors contain all qvector values which need to be calculated
        //////////////////////////////
        hid_t dspace_qv2 = H5Dget_space(ds_qv);

        hsize_t dfqt2[2];
        dfqt2[0]=finalqvectors.size(); dfqt2[1]=3;
        hid_t dspace_fqt2 = H5Screate_simple(2, dfqt2, NULL);             

        for(size_t i = 0; i < qindexes.size(); ++i)
        {
            // (over)write qvectors values for index positions
            hsize_t start[2];  // Start of hyperslab
            hsize_t count[2];  // Block count

            // source
            start[0]=i;start[1]=0;
            count[0]=1;count[1]=3;
            H5Sselect_hyperslab(dspace_fqt2,H5S_SELECT_SET,start,NULL,count,NULL);

            // target
            start[0]=qindexes[i];start[1]=0;
            count[0]=1;count[1]=3;
            H5Sselect_hyperslab(dspace_qv2,H5S_SELECT_SET,start,NULL,count,NULL);

            double* p_data = reinterpret_cast<double*>(const_cast<CartesianCoor3D*>(&finalqvectors[0]));
            H5Dwrite(ds_qv,H5T_NATIVE_DOUBLE,dspace_fqt2,dspace_qv2,H5P_DEFAULT,p_data);
        }
        H5Sclose(dspace_qv2);
        H5Sclose(dspace_fqt2);        
    }

    H5Dclose(ds_qv);
    H5Dclose(ds_fqt);
    H5Dclose(ds_fq);    
    H5Dclose(ds_checkpoint);
    
    H5Fclose(h5file);
    
    return qindexes;
}


std::vector<CartesianCoor3D> HDF5WriterService::get_qvectors(const std::vector<size_t>& qindexes) 
{    
    const std::string filename = m_filename;
    std::vector<CartesianCoor3D> qvectors;
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    
    hid_t ds_qv = H5Dopen(h5file,"qvectors",H5P_DEFAULT);
    hid_t ds_checkpoint = H5Dopen(h5file,"checkpoint",H5P_DEFAULT);
    
    hid_t dspace_qv = H5Dget_space(ds_qv);
    hid_t dspace_checkpoint = H5Dget_space(ds_checkpoint);
    
    hsize_t qvector_field[2];
    hsize_t qvector_field_max[2];
    H5Sget_simple_extent_dims(dspace_qv,qvector_field,qvector_field_max);

    hsize_t checkpoint_field[1];
    hsize_t checkpoint_field_max[1];
    H5Sget_simple_extent_dims(dspace_checkpoint,checkpoint_field,checkpoint_field_max);
    
    std::vector<CartesianCoor3D> h5qvectors(qvector_field[0]);
    H5Dread(ds_qv,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,reinterpret_cast<double*>(&h5qvectors[0]));
    
    for(size_t i = 0; i < qindexes.size(); ++i)
    {
        qvectors.push_back(h5qvectors[qindexes[i]]);
    }

    H5Sclose(dspace_qv);
    H5Sclose(dspace_checkpoint);

    H5Dclose(ds_qv);
    H5Dclose(ds_checkpoint);
    
    H5Fclose(h5file);

    return qvectors;
}


void HDF5WriterService::open() {
    m_hdf5_handle = H5Fopen(m_filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    m_ds_fqt = H5Dopen(m_hdf5_handle,"fqt",H5P_DEFAULT);
    m_ds_fq = H5Dopen(m_hdf5_handle,"fq",H5P_DEFAULT);    
    m_ds_checkpoint = H5Dopen(m_hdf5_handle,"checkpoint",H5P_DEFAULT);
}

void HDF5WriterService::close() {
    H5Dclose(m_ds_fqt);
    H5Dclose(m_ds_fq);    
    H5Dclose(m_ds_checkpoint);
    H5Fclose(m_hdf5_handle);
}

void HDF5WriterService::write(const  size_t qindex, const std::vector<complex<double> >& fqt,const std::complex<double>& fq)
{
    hid_t dspace_fqt = H5Dget_space(m_ds_fqt);
    hid_t dspace_fq = H5Dget_space(m_ds_fq);    
    hid_t dspace_checkpoint = H5Dget_space(m_ds_checkpoint);
    
    hsize_t fqt_start[3];  // Start of hyperslab
    hsize_t fqt_count[3];  // Block count
    fqt_start[0]=qindex;
    fqt_start[1]=0;
    fqt_start[2]=0;
    fqt_count[0]=1;
    fqt_count[1]=fqt.size();
    fqt_count[2]=2;

    H5Sselect_hyperslab(dspace_fqt,H5S_SELECT_SET,fqt_start,NULL,fqt_count,NULL);
    hsize_t dfqt1[3];
    dfqt1[0]=1;dfqt1[1]=fqt.size();dfqt1[2]=2; 
    hid_t dspace_fqt1 = H5Screate_simple(3, dfqt1, NULL);         
    double* p_data = reinterpret_cast<double*>(const_cast<complex<double>*>(&fqt[0]));
    H5Dwrite(m_ds_fqt,H5T_NATIVE_DOUBLE,dspace_fqt1,dspace_fqt,H5P_DEFAULT,p_data);

    H5Sclose(dspace_fqt1);

    // fq
    
    hsize_t fq_start[2];  // Start of hyperslab
    hsize_t fq_count[2];  // Block count
    fq_start[0]=qindex;
    fq_start[1]=0;
    fq_count[0]=1;
    fq_count[1]=2;

    H5Sselect_hyperslab(dspace_fq,H5S_SELECT_SET,fq_start,NULL,fq_count,NULL);
    hsize_t dfq1[2];
    dfq1[0]=1;dfq1[1]=2;
    hid_t dspace_fq1 = H5Screate_simple(2, dfq1, NULL);         
    double* p_data2 = reinterpret_cast<double*>(const_cast<complex<double>*>(&fq));
    H5Dwrite(m_ds_fq,H5T_NATIVE_DOUBLE,dspace_fq1,dspace_fq,H5P_DEFAULT,p_data2);

    H5Sclose(dspace_fq1);


    int ok =1 ;	    		                
                       
    hsize_t cp_start[1];  // Block count
    hsize_t cp_count[1];  // Block count
    cp_start[0]=qindex;
    cp_count[0]=1;
  
    H5Sselect_hyperslab(dspace_checkpoint,H5S_SELECT_SET,cp_start,NULL,cp_count,NULL);
    hsize_t dcp1[1];
    dcp1[0]=1; 
    hid_t dspace_cp1 = H5Screate_simple(1, dcp1, NULL);         
    H5Dwrite(m_ds_checkpoint,H5T_NATIVE_INT,dspace_cp1,dspace_checkpoint,H5P_DEFAULT,&ok);
    H5Sclose(dspace_cp1);
    
    
    H5Sclose(dspace_fqt);
    H5Sclose(dspace_fq);
    H5Sclose(dspace_checkpoint);
 
}


HDF5WriterService::HDF5WriterService(boost::asio::io_service& io_service,const string filename, const std::vector<CartesianCoor3D>& qvectors,size_t nf)
    : m_io_service(io_service),
      m_acceptor(io_service,boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(),0))
{
    m_filename = filename;
    m_listener = NULL;
    m_listener_status = false;
    // initialize and check file
    m_qindexes = init(qvectors,nf);
    m_qvectors = get_qvectors(m_qindexes);
    
    // spawn a dedicated thread for writing
    

}

void HDF5WriterService::run() {
    m_listener = new boost::thread(boost::bind(&HDF5WriterService::listener,this));
}

void HDF5WriterService::listener() {
    m_listener_status = true;
    
    boost::asio::ip::tcp::socket socket( m_io_service );

    while (true) {
        m_acceptor.accept(socket);

        HDF5WriterTag tag;
        size_t qindex; 
        size_t size;
        size_t count;        

        boost::asio::read(socket,boost::asio::buffer(&tag,sizeof(HDF5WriterTag))); 
        
        if (tag==HANGUP) {
            socket.close();
            break;
        }

        if (tag==WRITE) {

            boost::asio::read(socket,boost::asio::buffer(&count,sizeof(size_t))); 

            for(size_t i=0;i<count;i++) {
                boost::asio::read(socket,boost::asio::buffer(&qindex,sizeof(size_t))); 
                boost::asio::read(socket,boost::asio::buffer(&size,sizeof(size_t))); 

                std::vector<complex<double> >* p_data = new std::vector<complex<double> >(size/2);
                double* p_doubledata = (double*) &((*p_data)[0]);
                std::complex<double> fq;
                double* p_fq = (double*) &fq;  
                boost::asio::read(socket,boost::asio::buffer(p_doubledata,sizeof(double)*size)); 
                boost::asio::read(socket,boost::asio::buffer(p_fq,sizeof(double)*2)); 
  
                data_queue_fqt.push(make_pair(qindex,p_data));                
                data_queue_fq.push(make_pair(qindex,fq));
            }
        }        
        socket.close();
        
        if ((boost::posix_time::second_clock::universal_time()-m_lastflush) >
            (boost::posix_time::seconds(Params::Inst()->limits.times.iowrite_server)) ) {
            flush();
        }

        // use size from above. works since all elements have same size
        size_t data_bytesize = (sizeof(double)*2*size+sizeof(size_t)) * data_queue_fqt.size();
        data_bytesize += (sizeof(double)*2+sizeof(size_t)) * data_queue_fq.size();

        if (data_bytesize > Params::Inst()->limits.memory.iowrite_server) {
            flush();
        }
        
    }
    
    flush(); // make sure any remaining data in the buffer is written to disk
    
    m_listener_status = false;
}


void HDF5WriterService::flush() {
    m_lastflush = boost::posix_time::second_clock::universal_time();
    
    open();                

    while (data_queue_fqt.size()>0) {
        std::pair<size_t,std::vector<std::complex<double> >* > el = data_queue_fqt.front();
        data_queue_fqt.pop();
        size_t qindex = el.first;
        std::pair<size_t,std::complex<double> > el2 = data_queue_fq.front();
        data_queue_fq.pop();
        
        std::vector<std::complex<double> >* p_data = el.second;
        std::complex<double> fq = el2.second;
        write(qindex,*p_data,fq);        
        delete p_data;
    }

    close();
}


void HDF5WriterService::hangup() {
    
    if (m_listener!=NULL) {
        boost::asio::ip::tcp::socket socket( m_io_service );
        socket.connect(m_acceptor.local_endpoint());
        HDF5WriterTag tag = HANGUP;
        boost::asio::write(socket,boost::asio::buffer(&tag,sizeof(HDF5WriterTag))); 
        socket.close();
    }
    // block until listener routine has ended
    while (m_listener_status) {
        boost::this_thread::sleep(boost::posix_time::microseconds(100));
    }
    
}


void HDF5WriterClient::write(size_t qindex,const std::vector<std::complex<double> >& data) {
    
    if (!Params::Inst()->debug.iowrite.write) return;
    
    std::vector<std::complex<double> >* p_data = new std::vector<std::complex<double> >(data.size());
    *p_data = data;
    data_queue.push(make_pair(qindex,p_data));

    if ((boost::posix_time::second_clock::universal_time()-m_lastflush) >
        (boost::posix_time::seconds(Params::Inst()->limits.times.iowrite_client)) ){
        flush();
    }

    size_t data_bytesize = (sizeof(double)*2*data.size()+sizeof(size_t)) * data_queue.size();
    if (data_bytesize > Params::Inst()->limits.memory.iowrite_client) {
        flush();
    }
    
    if (!Params::Inst()->debug.iowrite.buffer) {
        flush();
    }
}


void HDF5WriterClient::flush() {

    m_lastflush = boost::posix_time::second_clock::universal_time();
    
    
    if (data_queue.size()>0) {
    
        // pre-compute fq
        std::queue<std::complex<double> > data_queue_fq;
        std::queue<std::pair<size_t,std::vector<std::complex<double> >* > > data_queue_fqt_copy;
        while (data_queue.size()>0) {
            std::pair<size_t,std::vector<std::complex<double> >* > el = data_queue.front();
            data_queue.pop();
            std::vector<std::complex<double> >* p_data = el.second;
            data_queue_fqt_copy.push(el);
            data_queue_fq.push(smath::reduce(*p_data));
        }

        boost::asio::io_service io_service;
        boost::asio::ip::tcp::socket socket( io_service );

        socket.connect(m_endpoint);
        
        HDF5WriterTag tag = WRITE;
        size_t count = data_queue_fqt_copy.size();
        
        boost::asio::write(socket,boost::asio::buffer(&tag,sizeof(HDF5WriterTag)));     
        boost::asio::write(socket,boost::asio::buffer(&count,sizeof(size_t)));     
        
        while (data_queue_fqt_copy.size()>0) {
            std::pair<size_t,std::vector<std::complex<double> >* > el = data_queue_fqt_copy.front();
            data_queue_fqt_copy.pop();
            complex<double> fq = data_queue_fq.front();
            data_queue_fq.pop();
            
            size_t qindex = el.first;
            std::vector<std::complex<double> >* p_data = el.second;

            size_t size = 2*p_data->size();
                        
            boost::asio::write(socket,boost::asio::buffer(&qindex,sizeof(size_t))); 
            boost::asio::write(socket,boost::asio::buffer(&size,sizeof(size_t))); 
            double* p_doubledata = (double*) &((*p_data)[0]);
            boost::asio::write(socket,boost::asio::buffer(p_doubledata,sizeof(double)*size)); 
            double* p_doubledata2 = (double*) &(fq);
            boost::asio::write(socket,boost::asio::buffer(p_doubledata2,sizeof(double)*2)); 
                        
            delete p_data;
        }
        
        socket.close();
    }
}

HDF5WriterClient::~HDF5WriterClient() {
    flush();
}
