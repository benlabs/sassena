/*
 *  This file is part of the software sassena
 *
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008-2010 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/data_stager.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/thread.hpp>

// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;

struct binary_max : public binary_function<size_t,size_t,size_t> {
    size_t operator() (size_t a,size_t b) {if (a>b) return a; else return b;}
};

DataStagerByFrame::DataStagerByFrame(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,DivAssignment assignment,Timer& timer) 
    : m_sample(sample),
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),    
    FC_assignment(assignment),
    timer_(timer)
{
    NN = allcomm_.size();  
    NNPP = partitioncomm_.size();
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->scattering.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = allcomm_.rank();

    // check data size requirements & allocate memory
    size_t data_bytesize = FC_assignment.size()*NA*3*sizeof(coor_t);
    size_t data_bytesize_indicator_max = 0;
    boost::mpi::all_reduce(partitioncomm_,data_bytesize,data_bytesize_indicator_max,binary_max());

    if (Params::Inst()->limits.stage.memory.data<data_bytesize_indicator_max) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
        }
        throw;
    }

    p_coordinates = (coor_t*) malloc(data_bytesize);
}

coor_t* DataStagerByFrame::stage() {
    
    timer_.start("st:first");
    stage_firstpartition();
    timer_.stop("st:first");

    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");
    
    timer_.start("st:fill");
    stage_fillpartitions();
    timer_.stop("st:fill");
    
    return p_coordinates;
}

void DataStagerByFrame::stage_firstpartition() {
    
    size_t LNF = FC_assignment.size();
    
    for(size_t f=0;f<LNF;f++) {
        timer_.start("st:load");
        CoordinateSet* p_cset = m_sample.coordinate_sets.load(FC_assignment[f]);
                    
        for(size_t n=0;n<NA;n++) {
            p_coordinates[f*3*NA+3*n]=p_cset->c1[n];
            p_coordinates[f*3*NA+3*n+1]=p_cset->c2[n];
            p_coordinates[f*3*NA+3*n+2]=p_cset->c3[n];            
        }
        timer_.stop("st:load");
    }
}

void DataStagerByFrame::stage_fillpartitions() {
    // create a communicator to broadcast between partitions.    
    boost::mpi::communicator interpartitioncomm_ = allcomm_.split(partitioncomm_.rank());
    boost::mpi::broadcast(interpartitioncomm_,p_coordinates,FC_assignment.size()*NA*3,0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// By Atom
////////////////////////////////////////////////////////////////////////////////////////////////////

DataStagerByAtom::DataStagerByAtom(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,DivAssignment assignment,Timer& timer) 
    : m_sample(sample),
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),    
    FC_assignment(assignment),
    timer_(timer)
{
    NN = allcomm_.size(); 
    NNPP = partitioncomm_.size();    
    NP = NN/NNPP;
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->scattering.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = allcomm_.rank();

    // check data size requirements & allocate memory
    size_t data_bytesize = FC_assignment.size()*NF*3*sizeof(coor_t);
    size_t data_bytesize_indicator_max = 0;
    boost::mpi::all_reduce(partitioncomm_,data_bytesize,data_bytesize_indicator_max,binary_max());

    if (Params::Inst()->limits.stage.memory.data<data_bytesize_indicator_max) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
        }
        throw;
    }

    p_coordinates = (coor_t*) malloc(data_bytesize);
}

coor_t* DataStagerByAtom::stage() {

    if (allcomm_.rank()==0) Info::Inst()->write("Staging first partition.");
    if (allcomm_.rank()<partitioncomm_.size()) {
        timer_.start("st:first");
        stage_firstpartition();
        timer_.stop("st:first");
    } 

    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");

    if (allcomm_.rank()==0) Info::Inst()->write("Staging remaining partitions.");

    timer_.start("st:fill");
    stage_fillpartitions();
    timer_.stop("st:fill");

    return p_coordinates;
    
//    stage_registration();    
//    stage_data();
//    return p_coordinates;
}

void DataStagerByAtom::stage_firstpartition() {
    size_t rank = partitioncomm_.rank();
    
    size_t buffer_bytesize = Params::Inst()->limits.stage.memory.buffer;
    size_t frame_bytesize = NA*3*sizeof(coor_t);
    size_t framesbuffer_maxsize = buffer_bytesize/frame_bytesize;
    
    if (framesbuffer_maxsize==0) {
        if (partitioncomm_.rank()==0) {
            Err::Inst()->write("Cannot load trajectory into buffer.");
            Err::Inst()->write(string("limits.memory.data_stager=")+boost::lexical_cast<string>(Params::Inst()->limits.stage.memory.buffer));
            Err::Inst()->write(string("requested=")+boost::lexical_cast<string>(frame_bytesize));            
        }
        throw;
    }
    
    if (partitioncomm_.rank()==0) {
        Info::Inst()->write(string("Initializing buffer size to: ")+boost::lexical_cast<string>(framesbuffer_maxsize));
    }
    coor_t* p_coordinates_buffer = (coor_t*) malloc(framesbuffer_maxsize*NA*3*sizeof(coor_t));
    std::vector< std::vector<size_t> > framesbuffer(NNPP);
    
    // align iteration with number of frames.
    size_t NFaligned = NF;
    if ((NF%NNPP)!=0) {
        size_t cycles = NF/NNPP;
        NFaligned = (cycles+1)*NNPP;
    }
    
    for(size_t f = 0; f < NFaligned; ++f)
    {
        size_t s;
        s = f%NNPP; // this is the responsible data server            
        
        if ( (rank==s) && (f<NF) ) {
            timer_.start("st:load");
            CoordinateSet* p_cset = m_sample.coordinate_sets.load(f);
            coor_t* p_data = &(p_coordinates_buffer[framesbuffer[s].size()*NA*3]);
            
            for(size_t n=0;n<NA;n++) {
                p_data[3*n]=p_cset->c1[n];
                p_data[3*n+1]=p_cset->c2[n];
                p_data[3*n+2]=p_cset->c3[n];            
            }
            delete p_cset;
            timer_.stop("st:load");
        }
        
        framesbuffer[s].push_back(f);
        
        if ( ((f+1)%NNPP) ==0 ) {
            if (framesbuffer[rank].size()==framesbuffer_maxsize) {
                timer_.start("st:distribute");
                distribute_coordinates(p_coordinates_buffer,framesbuffer,rank);
                timer_.stop("st:distribute");
                for(size_t i = 0; i < NNPP; ++i)
                {
                    framesbuffer[i].clear();                                                
                }
            }
        }
    }
    
    if (framesbuffer[rank].size()!=0) {
        timer_.start("st:distribute");
        distribute_coordinates(p_coordinates_buffer,framesbuffer,rank);
        timer_.stop("st:distribute");
        for(size_t i = 0; i < NNPP; ++i)
        {
            framesbuffer[i].clear();                                                
        }
    }
    
    timer_.start("st:wait");
    partitioncomm_.barrier();
    timer_.stop("st:wait");        
    
    free(p_coordinates_buffer);
}

void DataStagerByAtom::distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s) {
    
    size_t LNF = framesbuffer[s].size();

    // first nodes always gets maximum number of atoms in div assignment:
    DivAssignment zero_node_assignment(NNPP,0,NA);
    size_t maxatoms = zero_node_assignment.size();
    
    coor_t* p_alignedframe = (coor_t*) malloc(LNF*(maxatoms*NNPP)*3*sizeof(coor_t));

    // copy and align frames
    for(size_t i = 0; i < NNPP; ++i)
    {
        DivAssignment target_node_assignment(NNPP,i,NA);
        size_t off = target_node_assignment.offset();
        size_t len = target_node_assignment.size();
        
        for(size_t f = 0; f < LNF; ++f)
        {            
            coor_t* p_from = &(p_coordinates_buffer[f*NA*3+off*3]);
            coor_t* p_to = &(p_alignedframe[i*(maxatoms*LNF)*3 + i*maxatoms*3 ]);        
            memcpy(p_to,p_from,len*3*sizeof(coor_t));
        }
    }

    // exchange here
    coor_t* p_alignedframeOUT = (coor_t*) malloc(LNF*(maxatoms*NNPP)*3*sizeof(coor_t));
    boost::mpi::all_to_all(partitioncomm_,p_alignedframe,3*maxatoms*LNF,p_alignedframeOUT);
    free(p_alignedframe);

    // re-copy coordinates into right place
    size_t len = FC_assignment.size();
    for(size_t i = 0; i < NNPP; ++i)
    {
        for(size_t f = 0; f < LNF; ++f)
        {
            size_t firstframe = framesbuffer[0][f];
            size_t targetframe = firstframe+i;
            if (targetframe>NF) break;
            
            coor_t* p_from = &(p_alignedframeOUT[i*(maxatoms*LNF)*3]);
            for(size_t n = 0; n < len; ++n)
            {
                p_coordinates[ NF*n*3 + targetframe*3 ]     = p_from[ i*maxatoms*3 + n*3 ];
                p_coordinates[ NF*n*3 + targetframe*3 + 1 ] = p_from[ i*maxatoms*3 + n*3 + 1 ];
                p_coordinates[ NF*n*3 + targetframe*3 + 2 ] = p_from[ i*maxatoms*3 + n*3 + 2 ];
            }
        }
    }
    
    free(p_alignedframeOUT);    
}

void DataStagerByAtom::stage_fillpartitions() {    
    // create a communicator to broadcast between partitions.    
    boost::mpi::communicator interpartitioncomm_ = allcomm_.split(partitioncomm_.rank());
    boost::mpi::broadcast(interpartitioncomm_,p_coordinates,FC_assignment.size()*NF*3,0);
}

// end of file
