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
    H5::H5File h5file( filename,H5F_ACC_TRUNC);   

    H5::DataSpace fspace,mspace;
    int fill_val_int = 0;
    double fill_val_double = 0;
    
    hsize_t dims1[2],maxdims1[2],cdims1[2];
    dims1[0]=qvectors.size();
    dims1[1]=3;
    maxdims1[0]=H5S_UNLIMITED;
    maxdims1[1]=3;
    H5::DSetCreatPropList qvector_cparms;
    cdims1[0]=1;
    cdims1[1]=1;
    qvector_cparms.setChunk(2,cdims1);
    qvector_cparms.setFillValue( H5::PredType::NATIVE_DOUBLE, &fill_val_double);
    H5::DataSet ds_qv = h5file.createDataSet( "qvectors", H5::DataType(H5::PredType::NATIVE_DOUBLE), H5::DataSpace(2,dims1,maxdims1),qvector_cparms );

    hsize_t dims2[3],maxdims2[3],cdims2[3];
    dims2[0]=qvectors.size();
    dims2[1]=nf;
    dims2[2]=2;
    maxdims2[0]=H5S_UNLIMITED;
    maxdims2[1]=nf;
    maxdims2[2]=2;
    H5::DSetCreatPropList fqt_cparms;
    cdims2[0]=1;
    cdims2[1]=1;
    cdims2[2]=1;
    fqt_cparms.setChunk(3,cdims2);
    fqt_cparms.setFillValue( H5::PredType::NATIVE_DOUBLE, &fill_val_double);
    h5file.createDataSet( "fqt", H5::DataType(H5::PredType::NATIVE_DOUBLE), H5::DataSpace(3,dims2,maxdims2),fqt_cparms );
    
    hsize_t dims3[1],maxdims3[1],cdims3[1];
    dims3[0]=qvectors.size();
    maxdims3[0]=H5S_UNLIMITED;
    H5::DSetCreatPropList checkpoint_cparms;
    cdims3[0]=1;
    checkpoint_cparms.setChunk(1,cdims3);
    checkpoint_cparms.setFillValue( H5::PredType::NATIVE_INT, &fill_val_int);
    h5file.createDataSet( "checkpoint", H5::DataType(H5::PredType::NATIVE_INT), H5::DataSpace(1,dims3,maxdims3),checkpoint_cparms);
 
    std::vector<size_t> qindexes;   
    for(size_t i = 0; i < qvectors.size(); ++i) qindexes.push_back(i);

    // (over)write qvectors values for index positions
    hsize_t foffset1[2];  // Start of hyperslab
    foffset1[0]=0;
    foffset1[1]=0;
    fspace = ds_qv.getSpace();
    fspace.selectHyperslab(H5S_SELECT_SET,dims1,foffset1);
    mspace = H5::DataSpace(2,dims1);
    ds_qv.write(reinterpret_cast<double*>(const_cast<CartesianCoor3D*>(&qvectors[0])), H5::PredType::NATIVE_DOUBLE,mspace,fspace);
    h5file.close();
    
    return qindexes;
}


std::vector<size_t> init_reuse(const std::string filename,const std::vector<CartesianCoor3D>& qvectors,size_t nf)
{     
    //////////////////////////////
    // read info from data file
    //////////////////////////////
    H5::H5File h5file;
    h5file.openFile( filename,H5F_ACC_RDWR);
    
    H5::DataSet ds_qv = h5file.openDataSet("qvectors");
    H5::DataSet ds_fqt = h5file.openDataSet("fqt");
    H5::DataSet ds_checkpoint = h5file.openDataSet("checkpoint");
    
    hsize_t qvector_field[2];
    ds_qv.getSpace().getSimpleExtentDims(qvector_field);
    
    hsize_t checkpoint_field[1];
    ds_checkpoint.getSpace().getSimpleExtentDims(checkpoint_field);

    hsize_t fqt_field[3];
    ds_fqt.getSpace().getSimpleExtentDims(fqt_field);
    
    if (fqt_field[1]!=nf) {
        // this in an "incompatible" data file
        // mv old one to a backup
        int n=0;
        while (boost::filesystem::exists(filename+".backup-"+boost::lexical_cast<string>(n))) n++;
        boost::filesystem::rename(filename,filename+".backup"+boost::lexical_cast<string>(n));
        h5file.close();
        // escape by returning init_new
        return init_new(filename,qvectors,nf);
    }
    
    std::vector<CartesianCoor3D> h5qvectors(qvector_field[0]);
    ds_qv.read(reinterpret_cast<double*>(&h5qvectors[0]), H5::PredType::NATIVE_DOUBLE);
    
    std::vector<int> checkpoint(checkpoint_field[0]);
    ds_checkpoint.read(reinterpret_cast<int*>(&checkpoint[0]), H5::PredType::NATIVE_INT);
    
    //////////////////////////////
    // detect faulty entries
    // re-use those entries which have an ok of 0
    //////////////////////////////
    std::vector<size_t> qindexes;
    std::vector<CartesianCoor3D> oldqvectors;
    
    size_t faultyentries=0;
    for(size_t i = 0; i < checkpoint.size(); ++i)
    {
        if (checkpoint[i]==0) {
            qindexes.push_back(i);                 
            faultyentries++;
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
            if (fabs(qvectors[i].x-oldqvectors[j].x)>numeric_limits<double>::min()) continue;
            if (fabs(qvectors[i].y-oldqvectors[j].y)>numeric_limits<double>::min()) continue;
            if (fabs(qvectors[i].z-oldqvectors[j].z)>numeric_limits<double>::min()) continue;
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

    size_t extendsize=0;
    if (finalqvectors.size()<faultyentries) {
        qindexes.resize(finalqvectors.size());
    } else {
        extendsize = finalqvectors.size()-faultyentries;                
    }
    
    if (extendsize>0) {
        for(size_t i = 0; i < extendsize; ++i)
        {
            qindexes.push_back(qvector_field[0]+i);
        }
        qvector_field[0]+=extendsize;
        ds_qv.extend(qvector_field);
        
        checkpoint_field[0]+=extendsize;
        ds_checkpoint.extend(checkpoint_field);
    
        hsize_t fqt_field[3];
        ds_fqt.getSpace().getSimpleExtentDims(fqt_field);                
        fqt_field[0]+=extendsize;
        ds_fqt.extend(fqt_field);
    }
    
    //////////////////////////////
    // qindexes contain the indexes for the qvector storage points
    // finalqvectors contain all qvector values which need to be calculated
    //////////////////////////////
    
    for(size_t i = 0; i < qindexes.size(); ++i)
    {
        // (over)write qvectors values for index positions
        hsize_t start[3];  // Start of hyperslab
        hsize_t count[3];  // Block count
        start[0]=qindexes[i];start[1]=0;
        count[0]=1;count[1]=3;
        H5::DataSpace dspace_qv = ds_qv.getSpace();
        dspace_qv.selectHyperslab(H5S_SELECT_SET,count,start);
        hsize_t dims[2],maxdims[2];
        dims[0]=1;dims[1]=3;
        maxdims[0]=H5S_UNLIMITED;maxdims[1]=3;
        
        H5::DataSpace mspace_qv(2,dims,maxdims);
        ds_qv.write(reinterpret_cast<double*>(&finalqvectors[i]), H5::PredType::NATIVE_DOUBLE,mspace_qv,dspace_qv);
        
    }
    h5file.close();
    
    return qindexes;
}


std::vector<size_t> init(const string filename, const std::vector<CartesianCoor3D>& qvectors,size_t nf) {
if (boost::filesystem::exists(filename)) {
 if (H5::H5File::isHdf5(filename)) {
     return init_reuse(filename,qvectors,nf);
  } else {
      Err::Inst()->write(filename + string(" does not to be a HDF5 data file. aborting..."));
      throw;
  }
 } else {
     return init_new(filename,qvectors,nf);
 }
}
std::vector<CartesianCoor3D> get_qvectors(const std::string filename,const std::vector<size_t>& qindexes) 
{    
    std::vector<CartesianCoor3D> qvectors;
    H5::H5File h5file;

    h5file.openFile( filename,H5F_ACC_RDONLY);
    
    H5::DataSet ds_qv = h5file.openDataSet("qvectors");
    H5::DataSet ds_checkpoint = h5file.openDataSet("checkpoint");

    hsize_t qvector_field[2];
    ds_qv.getSpace().getSimpleExtentDims(qvector_field);

    hsize_t checkpoint_field[1];
    ds_checkpoint.getSpace().getSimpleExtentDims(checkpoint_field);

    std::vector<CartesianCoor3D> h5qvectors(qvector_field[0]);
    ds_qv.read(reinterpret_cast<double*>(&h5qvectors[0]), H5::PredType::NATIVE_DOUBLE);
    
    for(size_t i = 0; i < qindexes.size(); ++i)
    {
        qvectors.push_back(h5qvectors[qindexes[i]]);
    }
    
    h5file.close();
    return qvectors;
}

void store(const std::string filename,const  size_t qindex, const std::vector<complex<double> >& fqt)
{
    H5::H5File h5file;
    
    h5file.openFile(filename,H5F_ACC_RDWR);

    // write fqt data to hdf5file
    H5::DataSet ds_fqt = h5file.openDataSet("fqt");

    hsize_t fqt_start[3];  // Start of hyperslab
    hsize_t fqt_count[3];  // Block count
    fqt_start[0]=qindex;
    fqt_start[1]=0;
    fqt_start[2]=0;
    fqt_count[0]=1;
    fqt_count[1]=fqt.size();
    fqt_count[2]=2;

    H5::DataSpace dspace_fqt = ds_fqt.getSpace();
    dspace_fqt.selectHyperslab(H5S_SELECT_SET,fqt_count,fqt_start);

    hsize_t dims1[3],maxdims1[3];
    dims1[0]=1;
    dims1[1]=fqt.size();
    dims1[2]=2;
    maxdims1[0]=H5S_UNLIMITED;
    maxdims1[1]=H5S_UNLIMITED;
    maxdims1[2]=2;
    
    H5::DataSpace mspace_fqt = H5::DataSpace(3,dims1,maxdims1);
    
    ds_fqt.write(reinterpret_cast<double*>(const_cast<double*>(&fqt[0].real())), H5::PredType::NATIVE_DOUBLE,mspace_fqt,dspace_fqt);
    
    H5::DataSet ds_checkpoint = h5file.openDataSet("checkpoint");
    H5::DataSpace dspace_checkpoint = ds_checkpoint.getSpace();

    hsize_t ok_mdims[1] = {1};
    int ok =1 ;	    		                
                    
    hsize_t ok_foffset[1];  // Block count
    ok_foffset[0]=qindex;
    dspace_checkpoint.selectHyperslab(H5S_SELECT_SET,ok_mdims,ok_foffset);

    ds_checkpoint.write(&ok, H5::PredType::NATIVE_INT,H5::DataSpace (1,ok_mdims),dspace_checkpoint);

    h5file.close();
}

}
