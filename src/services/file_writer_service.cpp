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

void HDF5WriterService::init(size_t nf) {

    if (nf<1) {
        Err::Inst()->write("Number of frames must not be negative or zero!");    
        throw;
    }

    if (boost::filesystem::exists(m_filename)) {
        htri_t h5trialcode = H5Fis_hdf5(m_filename.c_str());  
        if (h5trialcode==0) {
            Err::Inst()->write(m_filename + string(" does not to be a HDF5 data file. aborting..."));
            throw;
        } else if (h5trialcode<0){
            Err::Inst()->write("Error when testing for HDF5 character of data file. aborting...");
            throw;         
        }
    } else {
        init_new(nf);
    }
}



void HDF5WriterService::init_new(size_t nf)
{     
    const string filename = m_filename;
    hid_t h5file = H5Fcreate( filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);   

    double fill_val_double = 0;

    // qvectors
    hsize_t dims[2];hsize_t mdims[2];hsize_t cdims[2];
    dims[0]=0;              dims[1]=3;
    mdims[0]=H5S_UNLIMITED; mdims[1]=3;
    cdims[0]=Params::Inst()->limits.signal.chunksize;         cdims[1]=3;
    hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk( dcpl, 2, cdims);
    H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
    if (Params::Inst()->limits.signal.alloc_early) {
        H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY);
    } else {
        H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
    }
    hid_t dspace = H5Screate_simple(2, dims, mdims); 
    hid_t ds = H5Dcreate(h5file, "qvectors", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
    H5Pclose(lcpl);
    H5Pclose(dcpl);
    H5Pclose(dapl);
    H5Sclose(dspace);
    H5Dclose(ds);

    if (Params::Inst()->scattering.signal.fqt) {
        hsize_t dims[3];hsize_t mdims[3];hsize_t cdims[3];
        dims[0]=0;              dims[1]=nf;     dims[2]=2;
        mdims[0]=H5S_UNLIMITED; mdims[1]=nf;     mdims[2]=2;
        size_t fqt_chunksize = Params::Inst()->limits.signal.chunksize;
        size_t fqt_chunksize2 = ((nf>0) && (nf<fqt_chunksize)) ? nf : fqt_chunksize;
        size_t fqt_chunksize1 = fqt_chunksize/fqt_chunksize2;
        cdims[0]=fqt_chunksize1; cdims[1]=fqt_chunksize2; cdims[2]=2;
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk( dcpl, 3, cdims);
        H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
        if (Params::Inst()->limits.signal.alloc_early) {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY);
        } else {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
        }
        hid_t dspace = H5Screate_simple(3, dims, mdims); 
        hid_t ds = H5Dcreate(h5file, "fqt", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
        H5Pclose(lcpl);
        H5Pclose(dcpl);
        H5Pclose(dapl);
        H5Sclose(dspace);
        H5Dclose(ds);
    }
    
    if (Params::Inst()->scattering.signal.fq0) {
        hsize_t dims[2];hsize_t mdims[2];hsize_t cdims[2];
        dims[0]=0;              dims[1]=2;     
        mdims[0]=H5S_UNLIMITED; mdims[1]=2;    
        cdims[0]=Params::Inst()->limits.signal.chunksize; cdims[1]=2; 
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk( dcpl, 2, cdims);
        H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
        if (Params::Inst()->limits.signal.alloc_early) {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY);
        } else {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
        }
        hid_t dspace = H5Screate_simple(2, dims, mdims); 
        hid_t ds = H5Dcreate(h5file, "fq0", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
        H5Pclose(lcpl);
        H5Pclose(dcpl);
        H5Pclose(dapl);
        H5Sclose(dspace);
        H5Dclose(ds);
    }
    
    if (Params::Inst()->scattering.signal.fq) {
        hsize_t dims[2];hsize_t mdims[2];hsize_t cdims[2];
        dims[0]=0;              dims[1]=2;     
        mdims[0]=H5S_UNLIMITED; mdims[1]=2;    
        cdims[0]=Params::Inst()->limits.signal.chunksize; cdims[1]=2; 
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk( dcpl, 2, cdims);
        H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
        if (Params::Inst()->limits.signal.alloc_early) {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY);
        } else {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
        }
        hid_t dspace = H5Screate_simple(2, dims, mdims); 
        hid_t ds = H5Dcreate(h5file, "fq", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
        H5Pclose(lcpl);
        H5Pclose(dcpl);
        H5Pclose(dapl);
        H5Sclose(dspace);
        H5Dclose(ds);
    }
    
    if (Params::Inst()->scattering.signal.fq2) {
        hsize_t dims[2];hsize_t mdims[2];hsize_t cdims[2];
        dims[0]=0;              dims[1]=2;     
        mdims[0]=H5S_UNLIMITED; mdims[1]=2;    
        cdims[0]=Params::Inst()->limits.signal.chunksize; cdims[1]=2; 
        hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dapl = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk( dcpl, 2, cdims);
        H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &fill_val_double);
        if (Params::Inst()->limits.signal.alloc_early) {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_EARLY);
        } else {
            H5Pset_alloc_time(dcpl,H5D_ALLOC_TIME_DEFAULT);
        }
        hid_t dspace = H5Screate_simple(2, dims, mdims); 
        hid_t ds = H5Dcreate(h5file, "fq2", H5T_NATIVE_DOUBLE, dspace,lcpl,dcpl,dapl);
        H5Pclose(lcpl);
        H5Pclose(dcpl);
        H5Pclose(dapl);
        H5Sclose(dspace);
        H5Dclose(ds);
    }
    
    H5Fclose(h5file);
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

std::vector<CartesianCoor3D> HDF5WriterService::get_qvectors()
{
    herr_t status ;
    const std::string filename = m_filename;
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

HDF5WriterService::HDF5WriterService(boost::asio::io_service& io_service,const string filename,size_t nf)
    : m_io_service(io_service),
      m_acceptor(io_service,boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(),0))
{
    m_filename = filename;
    m_listener = NULL;
    m_listener_status = false;
    // initialize and check file
    init(nf);
    
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
                CartesianCoor3D qvector;
                boost::asio::read(socket,boost::asio::buffer(&qvector,sizeof(CartesianCoor3D))); 

                data_qvectors.push_back(qvector);
                
                if (Params::Inst()->scattering.signal.fqt) {
                    boost::asio::read(socket,boost::asio::buffer(&size,sizeof(size_t))); 
                    std::vector<complex<double> >* p_data = new std::vector<complex<double> >(size/2);
                    double* p_doubledata = (double*) &((*p_data)[0]);
                    boost::asio::read(socket,boost::asio::buffer(p_doubledata,sizeof(double)*size));
                    data_fqt.push_back(p_data);
                }
                if (Params::Inst()->scattering.signal.fq) {
                    std::complex<double> fq;
                    double* p_fq = (double*) &fq;  
                    boost::asio::read(socket,boost::asio::buffer(p_fq,sizeof(double)*2)); 
                    data_fq.push_back(fq);
                }
            }
        }        
        socket.close();
        
        if ((boost::posix_time::second_clock::universal_time()-m_lastflush) >
            (boost::posix_time::seconds(Params::Inst()->limits.times.iowrite_server)) ) {
            flush();
        }

        // use size from above. works since all elements have same size
        size_t data_bytesize = (sizeof(double)*2*size) * data_fqt.size();
        data_bytesize += (sizeof(double)*3) * data_qvectors.size();
        data_bytesize += (sizeof(double)*2) * data_fq.size();

        if (data_bytesize > Params::Inst()->limits.memory.iowrite_server) {
            flush();
        }
        
    }
    
    flush(); // make sure any remaining data in the buffer is written to disk
    
    m_listener_status = false;
}


void HDF5WriterService::flush() {
    m_lastflush = boost::posix_time::second_clock::universal_time();
    
    cout << "opening h5file " <<  boost::posix_time::second_clock::universal_time() << endl;
    hid_t h5file = H5Fopen(m_filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    cout << "flushing qvectors " <<  boost::posix_time::second_clock::universal_time() << endl;
    if (data_qvectors.size()>0) {
        size_t nq = data_qvectors.size();

        cout << "find qv " <<  boost::posix_time::second_clock::universal_time() << endl;        

        if (H5LTfind_dataset(h5file,"qvectors")==1) {
            cout << "foudn qv " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            hsize_t dims[2];hsize_t start[2];hsize_t count[2];
            H5LTget_dataset_info(h5file,"qvectors",dims,NULL,NULL); 
            
            dims[0]+=nq;
            hid_t ds = H5Dopen(h5file,"qvectors",H5P_DEFAULT);
            H5Dset_extent(ds,dims);
            hid_t dspace = H5Dget_space(ds);

            hsize_t mdims[2];
            mdims[0]=nq; mdims[1]=3;
            hid_t mspace = H5Screate_simple(2, mdims, NULL);             

            start[0]=dims[0]-data_qvectors.size();start[1]=0;
            count[0]=nq;count[1]=3;
            H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
            start[0]=0;start[1]=0;
            count[0]=nq;count[1]=3;            
            H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

            double* p_data = reinterpret_cast<double*>(&data_qvectors[0]);
            cout << "write qv " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_data);            
            cout << "wrote qv " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
        }        
    }
    data_qvectors.clear();        
    
    cout << "flushing fqt " <<  boost::posix_time::second_clock::universal_time() << endl;
    if (data_fqt.size()>0) {
        size_t nq = data_fqt.size();
        size_t nf = data_fqt[0]->size();
        
        cout << "find fqt " <<  boost::posix_time::second_clock::universal_time() << endl;        
        
        if (H5LTfind_dataset(h5file,"fqt")==1) {
            cout << "found fqt " <<  boost::posix_time::second_clock::universal_time() << endl;        

            double* p_fqtdata = (double*) malloc(sizeof(double)*2*nq*nf);
            for(size_t i = 0; i < nq; ++i)
            {
                memcpy(&p_fqtdata[i*nf*2],&((*data_fqt[i])[0]),sizeof(double)*2*nf);
            }
            hsize_t dims[3];hsize_t start[3];hsize_t count[3];
            H5LTget_dataset_info(h5file,"fqt",dims,NULL,NULL); 
            
            dims[0]+=nq;
            hid_t ds = H5Dopen(h5file,"fqt",H5P_DEFAULT);
            H5Dset_extent(ds,dims);
            hid_t dspace = H5Dget_space(ds);

            hsize_t mdims[3];
            mdims[0]=nq; mdims[1]=nf;mdims[2]=2;
            hid_t mspace = H5Screate_simple(3, mdims, NULL);             

            start[0]=dims[0]-nq;start[1]=0;start[2]=0;
            count[0]=nq;count[1]=nf;count[2]=2;
            H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
            start[0]=0;start[1]=0;start[2]=0;
            count[0]=nq;count[1]=nf;count[2]=2;            
            H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

            cout << "write fqt " <<  boost::posix_time::second_clock::universal_time() << endl;        

            H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_fqtdata);  
            cout << "wrote fqt " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            free(p_fqtdata);          
        }  
               
        cout << "flushing fq0 " <<  boost::posix_time::second_clock::universal_time() << endl;            
        cout << "find fq0 " <<  boost::posix_time::second_clock::universal_time() << endl;        
        
        if (H5LTfind_dataset(h5file,"fq0")==1) {
            cout << "found fq0 " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            double* p_fq0data = (double*) malloc(sizeof(double)*2*nq);
            
            for(size_t i = 0; i < nq; ++i)
            {
                memcpy(&p_fq0data[i*2],&((*data_fqt[i])[0]),sizeof(double)*2);
            }
                        
            hsize_t dims[2];hsize_t start[2];hsize_t count[2];
            H5LTget_dataset_info(h5file,"fq0",dims,NULL,NULL); 
            
            dims[0]+=nq;
            hid_t ds = H5Dopen(h5file,"fq0",H5P_DEFAULT);
            H5Dset_extent(ds,dims);
            hid_t dspace = H5Dget_space(ds);

            hsize_t mdims[2];
            mdims[0]=nq; mdims[1]=2;
            hid_t mspace = H5Screate_simple(2, mdims, NULL);             

            start[0]=dims[0]-nq;start[1]=0;
            count[0]=nq;count[1]=2;
            H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
            start[0]=0;start[1]=0;
            count[0]=nq;count[1]=2;            
            H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

            cout << "write fq0 " <<  boost::posix_time::second_clock::universal_time() << endl;        

            H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_fq0data);
            cout << "wrote fq0 " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            free(p_fq0data);       
        }  
    }
    for(size_t i = 0; i < data_fqt.size(); ++i) delete data_fqt[i];
    data_fqt.clear();
    
    cout << "flushing fq " <<  boost::posix_time::second_clock::universal_time() << endl;
    if (data_fq.size()>0) {
        size_t nq = data_fq.size();

        cout << "find fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
        if (H5LTfind_dataset(h5file,"fq")==1) {
            cout << "found fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

            hsize_t dims[2];hsize_t start[2];hsize_t count[2];
            cout << "getinfo fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
            H5LTget_dataset_info(h5file,"fq",dims,NULL,NULL);             
            cout << "gotinfo fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            dims[0]+=nq;
            hid_t ds = H5Dopen(h5file,"fq",H5P_DEFAULT);
            cout << "setextent fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            H5Dset_extent(ds,dims);
            cout << "setextent fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

            cout << "getspace fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

            hid_t dspace = H5Dget_space(ds);
            cout << "gotspace fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

            hsize_t mdims[2];
            mdims[0]=nq; mdims[1]=2;
            cout << "createspace fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            hid_t mspace = H5Screate_simple(2, mdims, NULL);             
            cout << "createdspace fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

            start[0]=dims[0]-nq;start[1]=0;
            count[0]=nq;count[1]=2;
            H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
            start[0]=0;start[1]=0;
            count[0]=nq;count[1]=2;            
            H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);
            double* p_data = reinterpret_cast<double*>(&data_fq[0]);
            cout << "write fq " <<  boost::posix_time::second_clock::universal_time() << endl;        
            
            H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_data);
            cout << "wrote fq " <<  boost::posix_time::second_clock::universal_time() << endl;        

        }
        
        cout << "flushing fq2 " <<  boost::posix_time::second_clock::universal_time() << endl;
        if (H5LTfind_dataset(h5file,"fq2")==1) {
            hsize_t dims[2];hsize_t start[2];hsize_t count[2];
            H5LTget_dataset_info(h5file,"fq2",dims,NULL,NULL);             
            
            std::vector<std::complex<double> > fq2data(nq);
            
            for(size_t i = 0; i < nq; ++i)
            {
                fq2data[i]=data_fq[i]*conj(data_fq[i]);
            }
            
            dims[0]+=nq;
            hid_t ds = H5Dopen(h5file,"fq2",H5P_DEFAULT);
            H5Dset_extent(ds,dims);
            hid_t dspace = H5Dget_space(ds);

            hsize_t mdims[2];
            mdims[0]=nq; mdims[1]=2;
            hid_t mspace = H5Screate_simple(2, mdims, NULL);             

            start[0]=dims[0]-nq;start[1]=0;
            count[0]=nq;count[1]=2;
            H5Sselect_hyperslab(dspace,H5S_SELECT_SET,start,NULL,count,NULL);
            start[0]=0;start[1]=0;
            count[0]=nq;count[1]=2;            
            H5Sselect_hyperslab(mspace,H5S_SELECT_SET,start,NULL,count,NULL);

            double* p_data = reinterpret_cast<double*>(&fq2data[0]);
            H5Dwrite(ds,H5T_NATIVE_DOUBLE,mspace,dspace,H5P_DEFAULT,p_data);
            
        }   
    }    
    data_fq.clear();
    
    cout << "closing h5file " <<  boost::posix_time::second_clock::universal_time() << endl;
    H5Fclose(h5file);
    cout << "done h5file " <<  boost::posix_time::second_clock::universal_time() << endl;
    
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


void HDF5WriterClient::write(CartesianCoor3D qvector,const std::vector<std::complex<double> >& data) {
    
    if (!Params::Inst()->debug.iowrite.write) return;
    
    std::vector<std::complex<double> >* p_data = new std::vector<std::complex<double> >(data.size());
    *p_data = data;
    data_queue.push(make_pair(qvector,p_data));

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
        std::queue<std::pair<CartesianCoor3D,std::vector<std::complex<double> >* > > data_queue_fqt_copy;
        while (data_queue.size()>0) {
            std::pair<CartesianCoor3D,std::vector<std::complex<double> >* > el = data_queue.front();
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
            std::pair<CartesianCoor3D,std::vector<std::complex<double> >* > el = data_queue_fqt_copy.front();
            data_queue_fqt_copy.pop();
            complex<double> fq = data_queue_fq.front();
            data_queue_fq.pop();
            
            CartesianCoor3D qvector = el.first;
            std::vector<std::complex<double> >* p_data = el.second;

            size_t size = 2*p_data->size();
                        
            boost::asio::write(socket,boost::asio::buffer(&qvector,sizeof(CartesianCoor3D)));
            if (Params::Inst()->scattering.signal.fqt) {
                boost::asio::write(socket,boost::asio::buffer(&size,sizeof(size_t))); 
                double* p_doubledata = (double*) &((*p_data)[0]);
                boost::asio::write(socket,boost::asio::buffer(p_doubledata,sizeof(double)*size));                 
            } 
            if (Params::Inst()->scattering.signal.fq) {
                double* p_doubledata2 = (double*) &(fq);
                boost::asio::write(socket,boost::asio::buffer(p_doubledata2,sizeof(double)*2));                 
            }
                        
            delete p_data;
        }
        
        socket.close();
    }
}

HDF5WriterClient::~HDF5WriterClient() {
    flush();
}
