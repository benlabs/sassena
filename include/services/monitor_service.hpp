/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */
#ifndef IO__MONITOR_SERVICE_HPP_
#define IO__MONITOR_SERVICE_HPP_

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

#include <report/timer.hpp>


// other headers
#include "math/coor3d.hpp"

enum MonitorTag {MONITOR_HANGUP,MONITOR_UPDATE,MONITOR_RESET,MONITOR_SAMPLINGFACTOR};

class MonitorService {
    boost::asio::io_service& m_io_service;
    boost::asio::ip::tcp::acceptor m_acceptor;
    
    double m_from;
    double m_to;
    double m_current;
    
    Timer m_timer;
    
    std::string m_timerlabel;
    boost::thread* m_listener;
    bool m_listener_status;
    
    std::queue<double> update_thresholds;
    
    void listener();
    void update();
    void print();
    void reset_timer();
    void reset_state();

    size_t samplingfactor_;

public:
    MonitorService(boost::asio::io_service& io_service,double from,double to);
    
    boost::asio::ip::tcp::endpoint get_endpoint() { return m_acceptor.local_endpoint(); }
    
    void reset();
    
    void hangup();
    void run();
};

class MonitorClient {
    boost::asio::ip::tcp::endpoint m_endpoint;
    
    std::queue<double> update_thresholds;
    
    boost::posix_time::ptime lastupdate_;  
    size_t samplingfactor_;
public:
    MonitorClient(boost::asio::ip::tcp::endpoint server);
    void reset_server();
    
    void set_samplingfactor(size_t samplingfactor);
    void set_samplingfactor_server(size_t samplingfactor);
        
    void update(size_t rank,double progress);
};


#endif
