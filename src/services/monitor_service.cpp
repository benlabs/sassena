/** \file
This file implements a monitoring service, where clients push updates towards to head node, which then presents the user with up-to-date information about the progress.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "services/monitor_service.hpp"
#include <assert.h>

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/random/uniform_int.hpp>
#include <log.hpp>
#include <control.hpp>
#include <report/timer.hpp>

using namespace std;

MonitorService::MonitorService(boost::asio::io_service& io_service,double from,double to)
    : m_io_service(io_service),
      m_acceptor(io_service,boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(),0)),
      m_from(from),
      m_to(to)
{
   m_listener = NULL;
   m_listener_status = false;
   m_timerlabel = "progress";
   
   reset_state();
}

void MonitorService::run() {
    m_listener = new boost::thread(boost::bind(&MonitorService::listener,this));
    // triggers the internal timer
    m_timer.start(m_timerlabel);
}

void MonitorService::listener() {
    m_listener_status = true;
    
    boost::asio::ip::tcp::socket socket( m_io_service );
    
    std::map<size_t,double> progresses;
    
    while (true) {
        m_acceptor.accept(socket);

        update();
        MonitorTag tag;
        size_t rank; 
        double progress;

        socket.read_some(boost::asio::buffer(&tag,sizeof(MonitorTag)));
        
        if (tag==MONITOR_HANGUP) {
            break;
        }
        
        if (tag==MONITOR_RESET) {
            reset();
        }
        
        if (tag==MONITOR_SAMPLINGFACTOR) {
            socket.read_some(boost::asio::buffer(&samplingfactor_,sizeof(size_t)));            
        }
        
        if (tag==MONITOR_UPDATE) {
            socket.read_some(boost::asio::buffer(&rank,sizeof(size_t)));
            socket.read_some(boost::asio::buffer(&progress,sizeof(double)));
            
            if (progress>progresses[rank]) progresses[rank]=progress;

            double sum=0;
            for(std::map<size_t,double>::iterator pi = progresses.begin();pi!=progresses.end();pi++) {
                sum+=pi->second;
            }
            m_current = sum*samplingfactor_;
            update();
        }		
        
        socket.close();
    }

    update();
    m_listener_status = false;
}

void MonitorService::hangup() {
    if (m_listener!=NULL) {
        boost::asio::ip::tcp::socket socket( m_io_service );
        socket.connect(m_acceptor.local_endpoint());
        MonitorTag tag = MONITOR_HANGUP;
        socket.write_some(boost::asio::buffer(&tag,sizeof(MonitorTag))); 
        socket.close();
    }
    // block until listener routine has ended
    while (m_listener_status) {
        boost::this_thread::sleep(boost::posix_time::microseconds(100));
    }
}

void MonitorService::reset() {
    reset_state();
    reset_timer();
}

void MonitorService::reset_state() {
   double m_scale = m_to-m_from;
    m_current = m_from;
     // create a queue of update thresholds
    while (!update_thresholds.empty()) { update_thresholds.pop(); }
    update_thresholds.push(0.0001*m_scale);
    update_thresholds.push(0.001*m_scale);
    update_thresholds.push(0.005*m_scale);
    update_thresholds.push(0.01*m_scale);
    update_thresholds.push(0.02*m_scale);
    update_thresholds.push(0.03*m_scale);
    update_thresholds.push(0.04*m_scale);
    update_thresholds.push(0.05*m_scale);
    update_thresholds.push(0.1*m_scale);    
    update_thresholds.push(0.2*m_scale);
    update_thresholds.push(0.3*m_scale);
    update_thresholds.push(0.4*m_scale);
    update_thresholds.push(0.5*m_scale);
    update_thresholds.push(0.6*m_scale);
    update_thresholds.push(0.7*m_scale);
    update_thresholds.push(0.8*m_scale);
    update_thresholds.push(0.9*m_scale);
    update_thresholds.push(0.9*m_scale);
    update_thresholds.push(0.95*m_scale);
    update_thresholds.push(0.96*m_scale);
    update_thresholds.push(0.97*m_scale);
    update_thresholds.push(0.98*m_scale);
    update_thresholds.push(0.99*m_scale);
    update_thresholds.push(0.99*m_scale);
    update_thresholds.push(0.99*m_scale);
    update_thresholds.push(0.995*m_scale);
    update_thresholds.push(0.999*m_scale);
    update_thresholds.push(0.9999*m_scale);
    update_thresholds.push(1.0*m_scale);
}

void MonitorService::reset_timer() {
    m_timer.stop(m_timerlabel);
    m_timer.clear();
    m_timer.start(m_timerlabel);
}

void MonitorService::update() {
    
    m_timer.stop(m_timerlabel);
	m_timer.start(m_timerlabel);
	
    size_t oldsize = update_thresholds.size();

    if (oldsize<1) return;

    while(update_thresholds.front()<m_current) {
        update_thresholds.pop();
        if (update_thresholds.size()<1) break;
    }

    if (update_thresholds.size()!=oldsize) {
        print();
    }
}

void MonitorService::print() {
    double fraction_done = (m_current-m_from) / (m_to-m_from);
	//average time since start: 
	double totaltime = m_timer.sum(m_timerlabel);
	
	double etotal = ( 1.0/fraction_done) * totaltime;
	double eta = etotal - totaltime;

    stringstream ifraction_done_strstr;
    ifraction_done_strstr.precision(3);
    ifraction_done_strstr.setf(std::ios::fixed);
    ifraction_done_strstr << 100*fraction_done;
    
	Info::Inst()->write(string("Progress: ") + ifraction_done_strstr.str() 
	    + string("%, ETOTAL: ")+boost::lexical_cast<string>(etotal) 
	    + string(" s, ETA: ")+boost::lexical_cast<string>(eta)+ string(" s"));
    
}

MonitorClient::MonitorClient(boost::asio::ip::tcp::endpoint server) : m_endpoint(server) {
    // create a queue of update thresholds
    update_thresholds.push(0.0001);
    update_thresholds.push(0.001);
    update_thresholds.push(0.005);
    double p=0.01;
    while (p<1.0) {
        update_thresholds.push(p);        
        p+=0.01;
    }
    update_thresholds.push(1.0);
    
    lastupdate_ = boost::posix_time::second_clock::universal_time();
    samplingfactor_ = 1;
    
}

void MonitorClient::update(size_t rank,double progress) {
    
    if (!Params::Inst()->debug.monitor.update) return;
    
    // first test sampling criteria       
    if ((rank%samplingfactor_)!=0) return;
    
    size_t oldsize = update_thresholds.size();

    if (oldsize<1) return;

    while(update_thresholds.front()<progress) {
        update_thresholds.pop();
        if (update_thresholds.size()<1) break;
    }
        
    if (update_thresholds.size()!=0) {    
        // then test for minimum time delay
        if ((boost::posix_time::second_clock::universal_time()-lastupdate_) <=
            (boost::posix_time::seconds(Params::Inst()->limits.services.monitor.delay)) )
        {
                return;
        }
    }
        
    if (update_thresholds.size()!=oldsize) {
        // setup monitoring service
        boost::asio::io_service io_service;
        boost::asio::ip::tcp::socket socket( io_service );

        lastupdate_ = boost::posix_time::second_clock::universal_time();        
        try {
            socket.connect(m_endpoint);
            MonitorTag tag = MONITOR_UPDATE;
            socket.write_some(boost::asio::buffer(&tag,sizeof(MonitorTag)));        
            socket.write_some(boost::asio::buffer(&rank,sizeof(size_t)));
            socket.write_some(boost::asio::buffer(&progress,sizeof(double)));
            socket.close();            
        } catch(...) {
            Warn::Inst()->write("Unable to send update to monitor server");
            Warn::Inst()->write("Increase limits.services.monitor.delay and/or limits.services.monitor.sampling");
            Warn::Inst()->write(string("Current sampling factor: ")+boost::lexical_cast<string>(samplingfactor_));            
        }
    }
}

void MonitorClient::reset_server() {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket( io_service );
        
    socket.connect(m_endpoint);
    MonitorTag tag = MONITOR_RESET;
    socket.write_some(boost::asio::buffer(&tag,sizeof(MonitorTag)));        
    socket.close();
}

void MonitorClient::set_samplingfactor(size_t samplingfactor) {
    samplingfactor_=samplingfactor;
}

void MonitorClient::set_samplingfactor_server(size_t samplingfactor) {
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::socket socket( io_service );
        
    socket.connect(m_endpoint);
    MonitorTag tag = MONITOR_SAMPLINGFACTOR;
    socket.write_some(boost::asio::buffer(&tag,sizeof(MonitorTag)));        
    socket.write_some(boost::asio::buffer(&samplingfactor,sizeof(size_t)));            
    socket.close();
}



// end of file
