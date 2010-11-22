/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
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
    boost::mpi::communicator& allcomm_;
    boost::mpi::communicator& partitioncomm_;
    DivAssignment FC_assignment;
    
    Timer& timer_;
    
    size_t NFN;
    size_t NN;
    size_t NNPP;    
    size_t NP;    
    size_t NA;
    size_t NF;
        
    coor_t* p_coordinates;
//    void stage_registration();
//    void stage_data();
    void stage_firstpartition();
    void distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s);
    void stage_fillpartitions();
    
public:
    DataStagerByFrame(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm, DivAssignment assignment,Timer& timer);
    coor_t* stage();
    
};


class DataStagerByAtom  {
    Sample& m_sample;
    boost::mpi::communicator& allcomm_;
    boost::mpi::communicator& partitioncomm_;
    DivAssignment FC_assignment;
    
    Timer& timer_;
    
    size_t NFN;
    size_t NN;
    size_t NNPP;  
    size_t NP;  
    size_t NA;
    size_t NF;
    
    //owned by all clients, init in constructor
//    std::map<size_t, std::set<size_t> > FS_assignment_table;
//    std::map<size_t,std::vector<size_t> > FS_to_framelist_table;

//    std::vector<std::pair<std::pair<size_t,size_t>,std::vector<size_t> > > AtomOffLen_to_FClist;
    
    void distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s);
    void fill_coordinates(coor_t* p_localdata,size_t len,std::vector<size_t> frames);
//    void add_client(size_t off,size_t len, size_t client);
    
    coor_t* p_coordinates;
//    void stage_registration();
//    void stage_data();
    
    void stage_firstpartition();
    void stage_fillpartitions();
    
public:
    DataStagerByAtom(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm, DivAssignment assignment,Timer& timer);
  
    coor_t* stage();    
};

#endif