/*
 *  file_writer_service.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef IO__FILE_WRITER_SERVICE_HPP_
#define IO__FILE_WRITER_SERVICE_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <queue>
#include <vector>

// special library headers
#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <hdf5.h>


// other headers
#include "math/coor3d.hpp"

enum HDF5WriterTag {HANGUP,WRITE};

class HDF5WriterClient {
    boost::asio::ip::tcp::endpoint m_endpoint;
    boost::posix_time::ptime m_lastflush;
    
    std::queue<std::pair<size_t,std::vector<std::complex<double> >*> > data_queue;
    
public:
    HDF5WriterClient(boost::asio::ip::tcp::endpoint server) : m_endpoint(server) {}
    ~HDF5WriterClient();

    void write(size_t qindex,const std::vector<std::complex<double> >& data);
    void flush();    
};

class HDF5WriterService  {
    
    std::string m_filename;
    
    // these handles are kept when the file is opened:
    hid_t m_hdf5_handle;
    hid_t m_ds_checkpoint;
    hid_t m_ds_fqt;
    hid_t m_ds_fq;
    
    boost::posix_time::ptime m_lastflush;

    boost::asio::io_service& m_io_service;
    boost::asio::ip::tcp::acceptor m_acceptor;

    boost::thread* m_listener;
    bool m_listener_status;
    std::vector<size_t> m_qindexes;
    std::vector<CartesianCoor3D> m_qvectors;

    std::vector<size_t> init(const std::vector<CartesianCoor3D>& qvectors,size_t nf);
    std::vector<size_t> init_new(const std::vector<CartesianCoor3D>& qvectors,size_t nf);
    std::vector<size_t> init_reuse(const std::vector<CartesianCoor3D>& qvectors,size_t nf);

    void open();
    void close();
    void write(const size_t qindex, const std::vector<std::complex<double> >& fqt,const std::complex<double>& fq);

    bool test_fqt_dim(size_t nf);
    std::vector<CartesianCoor3D> get_qvectors(const std::vector<size_t>& qindexes);
    
    void listener();
    
    void flush();
    
    std::queue<std::pair<size_t,std::vector<std::complex<double> >*> > data_queue_fqt;
    std::queue<std::pair<size_t,std::complex<double> > > data_queue_fq;
public:
    
    HDF5WriterService(boost::asio::io_service& io_service, const std::string filename,const std::vector<CartesianCoor3D>& qvectors,size_t nf);

    boost::asio::ip::tcp::endpoint get_endpoint() { return m_acceptor.local_endpoint();}
    std::vector<size_t> get_qindexes() { return m_qindexes;}
    std::vector<CartesianCoor3D> get_qvectors() { return m_qvectors;}
    
    void hangup();
    void run();
};



#endif
