/*
 *  data_stager.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTER_DEVICES__DATA_STAGER_HPP_
#define SCATTER_DEVICES__DATA_STAGER_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <sys/time.h>
#include <vector>

// special library headers
#include <boost/asio.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpi.hpp>

// other headers
#include "sample.hpp"
#include "services.hpp"
#include "math/coor3d.hpp"
#include "decomposition/assignment.hpp"
#include "scatter_devices/scatter_factors.hpp"
#include "report/timer.hpp"

class DataStagerByFrame {
    Sample& m_sample;
    boost::mpi::communicator& m_comm;
    Assignment FC_assignment;
    
    size_t NFN;
    size_t NN;
    size_t NA;
    size_t NF;
    
    //owned by all clients, init in constructor
    std::map<size_t, std::set<size_t> > FS_assignment_table;
    std::map<size_t,std::vector<size_t> > FS_to_framelist_table;

    // used by file servers
    std::map<size_t,std::vector<size_t> > Frame_to_FClist; 
    
    coor_t* p_coordinates;
    void stage_registration();
    void stage_data();

public:
    DataStagerByFrame(Sample& sample,boost::mpi::communicator& comm, Assignment assignment);
    coor_t* stage();
    
};


class DataStagerByAtom  {
    Sample& m_sample;
    boost::mpi::communicator& m_comm;
    Assignment FC_assignment;
    
    size_t NFN;
    size_t NN;
    size_t NA;
    size_t NF;
    
    //owned by all clients, init in constructor
    std::map<size_t, std::set<size_t> > FS_assignment_table;
    std::map<size_t,std::vector<size_t> > FS_to_framelist_table;

    std::vector<std::pair<std::pair<size_t,size_t>,std::vector<size_t> > > AtomOffLen_to_FClist;
    
    void distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s);
    void fill_coordinates(coor_t* p_localdata,size_t len,std::vector<size_t> frames);
    void add_client(size_t off,size_t len, size_t client);
    
    coor_t* p_coordinates;
    void stage_registration();
    void stage_data();
    
public:
    DataStagerByAtom(Sample& sample,boost::mpi::communicator& comm, Assignment assignment);
    
    coor_t* stage();    
};

#endif