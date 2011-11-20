/** \file 
This file contains a class which implements the data staging logic for the trajectory data. The first partition reads the trajectory data from disc and distributes the trajectory to any other partition. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/


#ifndef STAGER__DATA_STAGER_HPP_
#define STAGER__DATA_STAGER_HPP_

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
#include "report/timer.hpp"
#include "stager/coordinate_writer.hpp"

/** 
Loads coordinate data, performs frame decomposition of the trajectory data and places it efficiently into the distributed memory of the parallel partition by using div logic. 
*/

class DataStagerByFrame {
    Sample& m_sample;
    boost::mpi::communicator& allcomm_;
    boost::mpi::communicator& partitioncomm_;
    
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

    void write(std::string filename,std::string format);    
public:
    DataStagerByFrame(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm, Timer& timer);
    coor_t* stage();
    
};

/** 
Loads coordinate data, performs atom decomposition of the trajectory data and places it efficiently into the distributed memory of the parallel partition by using div logic. 
*/
class DataStagerByAtom  {
    Sample& m_sample;
    boost::mpi::communicator& allcomm_;
    boost::mpi::communicator& partitioncomm_;
    
    Timer& timer_;
    
    size_t NFN;
    size_t NN;
    size_t NNPP;  
    size_t NP;  
    size_t NA;
    size_t NF;
    
    void fill_coordinates(coor_t* p_localdata,size_t len,std::vector<size_t> frames);
    
    coor_t* p_coordinates;
    
    void stage_firstpartition();
    void stage_fillpartitions();

    void write(std::string filename,std::string format);
    
public:
    DataStagerByAtom(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,Timer& timer);
  
    coor_t* stage();    
};

#endif