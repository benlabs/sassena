/*
 *  h5_fqt_interface.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/io/h5_fqt_interface.hpp"

// standard header
#include <complex>
#include <fstream>
#include <limits>

// special library headers
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// other headers
#include "log.hpp"
#include "math/coor3d.hpp"


using namespace std;

namespace H5FQTInterface {

std::vector<size_t> init_new(const string filename,const std::vector<CartesianCoor3D>& qvectors,size_t nf)
{     
    hid_t h5file = H5Fcreate( filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);   

    int fill_val_int = 0;
    double fill_val_double = 0;
    
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
    H5Pset_alloc_time(qvector_dcpl,H5D_ALLOC_TIME_EARLY);
    hid_t dspace_qv = H5Screate_simple(2, dims1, maxdims1); 
    hid_t ds_qv = H5Dcreate(h5file, "qvectors", H5T_NATIVE_DOUBLE, dspace_qv, qvector_lcpl,qvector_dcpl,qvector_dapl);
    H5Pclose(qvector_lcpl);
    H5Pclose(qvector_dcpl);
    H5Pclose(qvector_dapl);
    H5Sclose(dspace_qv);
    H5Dclose(ds_qv);

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
    H5Pset_alloc_time(fqt_dcpl,H5D_ALLOC_TIME_EARLY);
    hid_t dspace_fqt = H5Screate_simple(3, dims2, maxdims2); 
    hid_t ds_fqt = H5Dcreate(h5file, "fqt", H5T_NATIVE_DOUBLE, dspace_fqt, fqt_lcpl,fqt_dcpl,fqt_dapl);
    H5Pclose(fqt_lcpl);
    H5Pclose(fqt_dcpl);
    H5Pclose(fqt_dapl);
    H5Sclose(dspace_fqt);
    H5Dclose(ds_fqt);
    
    hsize_t dims3[1],maxdims3[1],cdims3[1];
    dims3[0]=qvectors.size();
    maxdims3[0]=H5S_UNLIMITED;
    cdims3[0]=10000;
    hid_t cp_lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t cp_dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t cp_dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( cp_dcpl, 1, cdims3);
    H5Pset_fill_value( cp_dcpl, H5T_NATIVE_INT, &fill_val_int);
    H5Pset_alloc_time(cp_dcpl,H5D_ALLOC_TIME_EARLY);
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

bool test_fqt_dim(const std::string filename, size_t nf) {
    
    bool retval = false;
    //////////////////////////////
    // read info from data file
    //////////////////////////////
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
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


std::vector<size_t> init_reuse(const std::string filename,const std::vector<CartesianCoor3D>& qvectors,size_t nf)
{     
    if (test_fqt_dim(filename,nf)==false) {
        // this in an "incompatible" data file
        // mv old one to a backup
        int n=0;
        while (boost::filesystem::exists(filename+".backup-"+boost::lexical_cast<string>(n))) n++;
        boost::filesystem::rename(filename,filename+".backup-"+boost::lexical_cast<string>(n));
        // escape by returning init_new
        return init_new(filename,qvectors,nf);  
    };
    
    //////////////////////////////
    // read info from data file
    //////////////////////////////
    hid_t h5file = H5Fopen( filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t ds_qv = H5Dopen(h5file,"qvectors",H5P_DEFAULT);
    hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);
    hid_t ds_checkpoint = H5Dopen(h5file,"checkpoint",H5P_DEFAULT);
        
    hid_t dspace_qv = H5Dget_space(ds_qv);
    hid_t dspace_fqt = H5Dget_space(ds_fqt);
    hid_t dspace_checkpoint = H5Dget_space(ds_checkpoint);
    
    hsize_t qvector_field[2];
    hsize_t checkpoint_field[1];
    hsize_t fqt_field[3];
    hsize_t qvector_field_max[2];
    hsize_t checkpoint_field_max[1];
    hsize_t fqt_field_max[3];
    
    H5Sget_simple_extent_dims(dspace_qv,qvector_field,qvector_field_max);
    H5Sget_simple_extent_dims(dspace_fqt,fqt_field,fqt_field_max);
    H5Sget_simple_extent_dims(dspace_checkpoint,checkpoint_field,checkpoint_field_max);
    
    H5Sclose(dspace_qv);
    H5Sclose(dspace_fqt);
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
        }    

        //////////////////////////////
        // qindexes contain the indexes for the qvector storage points
        // finalqvectors contain all qvector values which need to be calculated
        //////////////////////////////
        hid_t dspace_qv2 = H5Dget_space(ds_qv);

        hsize_t dfq1[2];
        dfq1[0]=finalqvectors.size(); dfq1[1]=3;
        hid_t dspace_fq = H5Screate_simple(2, dfq1, NULL);             

        for(size_t i = 0; i < qindexes.size(); ++i)
        {
            // (over)write qvectors values for index positions
            hsize_t start[2];  // Start of hyperslab
            hsize_t count[2];  // Block count

            // source
            start[0]=i;start[1]=0;
            count[0]=1;count[1]=3;
            H5Sselect_hyperslab(dspace_fq,H5S_SELECT_SET,start,NULL,count,NULL);

            // target
            start[0]=qindexes[i];start[1]=0;
            count[0]=1;count[1]=3;
            H5Sselect_hyperslab(dspace_qv2,H5S_SELECT_SET,start,NULL,count,NULL);

            double* p_data = reinterpret_cast<double*>(const_cast<CartesianCoor3D*>(&finalqvectors[0]));
            H5Dwrite(ds_qv,H5T_NATIVE_DOUBLE,dspace_fq,dspace_qv2,H5P_DEFAULT,p_data);
        }
        H5Sclose(dspace_qv2);
        H5Sclose(dspace_fq);        
    }

    H5Dclose(ds_qv);
    H5Dclose(ds_fqt);
    H5Dclose(ds_checkpoint);
    
    H5Fclose(h5file);
    
    return qindexes;
}


std::vector<size_t> init(const string filename, const std::vector<CartesianCoor3D>& qvectors,size_t nf) {

if (nf<1) {
    Err::Inst()->write("Number of frames must not be negative or zero!");    
    throw;
}
    
    
if (boost::filesystem::exists(filename)) {
    htri_t h5trialcode = H5Fis_hdf5(filename.c_str());  
    if (h5trialcode>0) {
        return init_reuse(filename,qvectors,nf);
     } else if (h5trialcode==0) {
         Err::Inst()->write(filename + string(" does not to be a HDF5 data file. aborting..."));
         throw;
     } else {
         Err::Inst()->write("Error when testing for HDF5 character of data file. aborting...");
         throw;         
     }
    } else {
        return init_new(filename,qvectors,nf);
    }
}
std::vector<CartesianCoor3D> get_qvectors(const std::string filename,const std::vector<size_t>& qindexes) 
{    
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

void store(const std::string filename,const  size_t qindex, const std::vector<complex<double> >& fqt)
{
    hid_t h5file = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

    // write fqt data to hdf5file
    hid_t ds_fqt = H5Dopen(h5file,"fqt",H5P_DEFAULT);
    hid_t ds_checkpoint = H5Dopen(h5file,"checkpoint",H5P_DEFAULT);

    hid_t dspace_fqt = H5Dget_space(ds_fqt);
    hid_t dspace_checkpoint = H5Dget_space(ds_checkpoint);
    
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
    H5Dwrite(ds_fqt,H5T_NATIVE_DOUBLE,dspace_fqt1,dspace_fqt,H5P_DEFAULT,p_data);

    H5Sclose(dspace_fqt1);

    int ok =1 ;	    		                
                    
    hsize_t cp_start[1];  // Block count
    hsize_t cp_count[1];  // Block count
    cp_start[0]=qindex;
    cp_count[0]=1;
  
    H5Sselect_hyperslab(dspace_checkpoint,H5S_SELECT_SET,cp_start,NULL,cp_count,NULL);
    hsize_t dcp1[1];
    dcp1[0]=1; 
    hid_t dspace_cp1 = H5Screate_simple(1, dcp1, NULL);         
    H5Dwrite(ds_checkpoint,H5T_NATIVE_INT,dspace_cp1,dspace_checkpoint,H5P_DEFAULT,&ok);
    H5Sclose(dspace_cp1);
    
    
    H5Sclose(dspace_fqt);
    H5Sclose(dspace_checkpoint);

    H5Dclose(ds_fqt);
    H5Dclose(ds_checkpoint);

    H5Fclose(h5file);
    
}

}
