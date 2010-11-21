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
        
    // determine number of nodes which act as file servers:
    NFN = Params::Inst()->limits.stage.nodes;
    
    if (NFN>NN) NFN=NN; 
    if (NFN>NF) {
        NFN=NF;  
        if (allcomm_.rank()==0) {
            Info::Inst()->write("Number of data servers limited by number of frames");
        }
    }    
    
    if (NFN==0) {
        if (allcomm_.rank()==0) Info::Inst()->write(string("Setting limits.stage.nodes to partition size (0=automatic,limits.stage.nodes)"));
        if (allcomm_.rank()==0) Info::Inst()->write(string("limits.stage.nodes=")+boost::lexical_cast<string>(partitioncomm_.size()));
        NFN=partitioncomm_.size();    
    }
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
    size_t rank = allcomm_.rank();
    
    size_t buffer_bytesize = Params::Inst()->limits.stage.memory.buffer;
    size_t frame_bytesize = NA*3*sizeof(coor_t);
    size_t framesbuffer_maxsize = buffer_bytesize/frame_bytesize;
    
    if (framesbuffer_maxsize==0) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("Cannot load trajectory into buffer.");
            Err::Inst()->write(string("limits.stage.memory.buffer=")+boost::lexical_cast<string>(Params::Inst()->limits.stage.memory.buffer));
            Err::Inst()->write(string("requested=")+boost::lexical_cast<string>(frame_bytesize));            
        }
        throw;
    }

    size_t mode; // 0 = mod
    size_t modblock = Params::Inst()->limits.stage.modblock;
    if (Params::Inst()->limits.stage.mode=="mod") {
        if (allcomm_.rank()==0) Info::Inst()->write(string("limits.stage.mode=mod"));
        if (allcomm_.rank()==0) Info::Inst()->write(string("limits.stage.modblock=")+boost::lexical_cast<string>(modblock));                
        mode = 0;
    } else if (Params::Inst()->limits.stage.mode=="div") {
        if (allcomm_.rank()==0) {
            Info::Inst()->write(string("limits.stage.mode=div"));
            Info::Inst()->write(string("Adjusting parameters due to div logic:"));
            Info::Inst()->write(string("Setting limits.stage.memory.buffer to 1 frame"));
            Info::Inst()->write(string("Setting limits.stage.nodes=")+boost::lexical_cast<string>(partitioncomm_.size()));            
            Info::Inst()->write(string("Setting limits.stage.barrier=2"));            
        }
        framesbuffer_maxsize = 1;
        NFN = partitioncomm_.size();
        Params::Inst()->limits.stage.barrier = 2;
        mode = 1;
    } else {
        Err::Inst()->write(string("Stage mode not understood: ")+Params::Inst()->limits.stage.mode);
        Err::Inst()->write(string("limits.stage.mode= mod, div"));
        throw;
    }
    
    if (allcomm_.rank()==0) {
        Info::Inst()->write(string("Initializing buffer size to: ")+boost::lexical_cast<string>(framesbuffer_maxsize));
    }
    coor_t* p_coordinates_buffer = (coor_t*) malloc(framesbuffer_maxsize*NA*3*sizeof(coor_t));
    std::vector< std::vector<size_t> > framesbuffer(NFN);
    
    
    for(size_t f=0;f<NF;f++) {
        size_t s;
        if (mode==0) {
            s = (f/modblock)%NFN; // this is the responsible data server            
        } else {
            s = (f*NFN)/NF; // this is the responsible data server                        
        }
        
        if (framesbuffer[s].size()==framesbuffer_maxsize) {
            timer_.start("st:distribute");
            distribute_coordinates(p_coordinates_buffer,framesbuffer,s);     
            
            if (mode==0) {
                timer_.start("st:wait");
                allcomm_.barrier();
                timer_.stop("st:wait");                
            }
                   
            timer_.stop("st:distribute");
            framesbuffer[s].clear();
        }
        
        if (rank==s) {
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
    }
    
    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");
    
    for(size_t i = 0; i < framesbuffer.size(); ++i)
    {
        if (framesbuffer[i].size()!=0) {
            timer_.start("st:distribute");
            distribute_coordinates(p_coordinates_buffer,framesbuffer,i);      
            
            if (mode==0) {
                timer_.start("st:wait");
                allcomm_.barrier();
                timer_.stop("st:wait");                
            }     
            timer_.stop("st:distribute");
            
            framesbuffer[i].clear();
        }
    }
    
    free(p_coordinates_buffer);
}


void DataStagerByFrame::distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s) {
    
    size_t LNF = framesbuffer[s].size();
    size_t rank = allcomm_.rank();  
    
    size_t barrier = Params::Inst()->limits.stage.barrier;

    for(size_t i = 0; i < framesbuffer[s].size(); ++i)
    {
        size_t f = framesbuffer[s][i];
        size_t target_node = (f*NNPP)/NF;

        if (rank==s) {
            // send data
            if (target_node==s) {
                coor_t* p_localdata = &(p_coordinates[FC_assignment.index(f)*NA*3]);
                memcpy(p_localdata,p_coordinates_buffer,3*NA*sizeof(coor_t));                
            } else {
                allcomm_.send(target_node,0,p_coordinates_buffer,3*NA);
            }
        } else {
            if (rank==target_node) {
                coor_t* p_localdata = &(p_coordinates[FC_assignment.index(f)*NA*3]);
                allcomm_.recv(s,0,p_localdata,3*NA); 
            }
        }
       
        if ( ((i+1)%barrier) == 0) {
            timer_.start("st:wait");
            allcomm_.barrier();
            timer_.stop("st:wait");
        }
    }
}


void DataStagerByFrame::stage_fillpartitions() {
    
    bool firstpartition=true;
    if (allcomm_.rank()>=partitioncomm_.size()) firstpartition=false;
    size_t partitions = allcomm_.size()/partitioncomm_.size();

    if (firstpartition) {
        for(size_t i = 1; i < partitions; ++i)
        {
            size_t target_node = i*partitioncomm_.size() + partitioncomm_.rank();
            allcomm_.send(target_node,0,p_coordinates,FC_assignment.size()*NA*3);            
        }
    } else {
        size_t source_node = partitioncomm_.rank();
        allcomm_.recv(source_node,0,p_coordinates,FC_assignment.size()*NA*3);            
    }
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

//    Info::Inst()->write(string("Maximum memory allocated for coordinates (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
    p_coordinates = (coor_t*) malloc(data_bytesize);
        
    // determine number of nodes which act as file servers:
    NFN = Params::Inst()->limits.stage.nodes;
        
    if (NFN>NN) NFN=NN; 
    if (NFN>NF) {
        NFN=NF;  
        if (allcomm_.rank()==0) {
            Info::Inst()->write("Number of data servers limited by number of frames.");
        }
    }
    
    if (NFN==0) {
        if (allcomm_.rank()==0) Info::Inst()->write(string("Setting limits.stage.nodes to partition size (0=automatic,limits.stage.nodes)"));
        if (allcomm_.rank()==0) Info::Inst()->write(string("limits.stage.nodes=")+boost::lexical_cast<string>(partitioncomm_.size()));
        NFN=partitioncomm_.size();    
    }
}

coor_t* DataStagerByAtom::stage() {

    if (allcomm_.rank()==0) Info::Inst()->write("Staging first partition.");
    timer_.start("st:first");
    stage_firstpartition();
    timer_.stop("st:first");

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
    size_t rank = allcomm_.rank();
    
    // used by server
    //size_t frame_bytesize = NA*3*sizeof(double);
//    coor_t* p_coordinates_buffer = (coor_t*) malloc(NA*3*sizeof(coor_t));

    //std::set<size_t>& assigned_frames = FS_assignment_table[rank];

    // initialize local coordinates buffer
    
    size_t buffer_bytesize = Params::Inst()->limits.stage.memory.buffer;
    size_t frame_bytesize = NA*3*sizeof(coor_t);
    size_t framesbuffer_maxsize = buffer_bytesize/frame_bytesize;
    
    if (framesbuffer_maxsize==0) {
        if (allcomm_.rank()==0) {
            Err::Inst()->write("Cannot load trajectory into buffer.");
            Err::Inst()->write(string("limits.memory.data_stager=")+boost::lexical_cast<string>(Params::Inst()->limits.stage.memory.buffer));
            Err::Inst()->write(string("requested=")+boost::lexical_cast<string>(frame_bytesize));            
        }
        throw;
    }
    
    if (allcomm_.rank()==0) {
        Info::Inst()->write(string("Initializing buffer size to: ")+boost::lexical_cast<string>(framesbuffer_maxsize));
    }
    coor_t* p_coordinates_buffer = (coor_t*) malloc(framesbuffer_maxsize*NA*3*sizeof(coor_t));
    std::vector< std::vector<size_t> > framesbuffer(NFN);
    
    size_t modblock = Params::Inst()->limits.stage.modblock;
    if (Params::Inst()->limits.stage.mode=="mod") {
        if (allcomm_.rank()==0) {
            Info::Inst()->write(string("limits.stage.mode=mod"));
            Info::Inst()->write(string("limits.stage.modblock=")+boost::lexical_cast<string>(modblock));        
        }
    } else if (Params::Inst()->limits.stage.mode=="div") {
        Info::Inst()->write(string("limits.stage.mode=div"));
        Warn::Inst()->write(string("Staging in div logic not supporting for incoherent scattering"));
        Warn::Inst()->write(string("Reverting to mod logic: limits.stage.mode=mod"));
    } else {
        Err::Inst()->write(string("Stage mode not understood: ")+Params::Inst()->limits.stage.mode);
        Err::Inst()->write(string("limits.stage.mode= mod, div"));
        throw;
    }
    
    // align iteration with number of frames.
    size_t NFaligned = NF;
    if ((NF%NNPP)!=0) {
        size_t cycles = NF/NNPP;
        NFaligned += ( NF - cycles*NNPP);
    }
    
    for(size_t f = 0; f < NFaligned; ++f)
    {
        size_t s;
        s = (f/modblock)%NFN; // this is the responsible data server            
        
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
        
        // push frame on buffer
        framesbuffer[s].push_back(f);
        
        if ( ((NFaligned+1)%NNPP) ==0 ) {
            timer_.start("st:distribute");
            distribute_coordinates(p_coordinates_buffer,framesbuffer,rank);
            timer_.stop("st:distribute");
            framesbuffer.clear();            
        }
    }
    
    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");        
    
    free(p_coordinates_buffer);
}

void DataStagerByAtom::copyalign_frame(coor_t* p_alignedframe, coor_t* p_frame,size_t block) {
    for(size_t i = 0; i < NNPP; ++i)
    {
        DivAssignment target_node_assignment(NNPP,i,NA);
        size_t off = target_node_assignment.offset();
        size_t len = target_node_assignment.size();

        coor_t* p_from = &(p_frame[off*3]);        
        coor_t* p_to = &(p_alignedframe[i*block*3]);
        memcpy(p_to,p_from,len*3*sizeof(coor_t));
    }
}

void DataStagerByAtom::fill_alignedframe(coor_t* p_alignedatoms,size_t block,size_t firstframe) {
    
    size_t off = FC_assignment.offset();
    size_t len = FC_assignment.size();
    
    for(size_t i = 0; i < NNPP; ++i)
    {
        size_t targetframe = firstframe+i;
        if (targetframe>NF) break;
        
        for(size_t n = 0; n < len; ++n)
        {
            p_coordinates[ NF*n*3 + targetframe*3 ]     = p_alignedatoms[ i*block*3 + n*3 ];
            p_coordinates[ NF*n*3 + targetframe*3 + 1 ] = p_alignedatoms[ i*block*3 + n*3 + 1 ];
            p_coordinates[ NF*n*3 + targetframe*3 + 2 ] = p_alignedatoms[ i*block*3 + n*3 + 2 ];
        }
    }
}

void DataStagerByAtom::distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s) {
    
    size_t LNF = framesbuffer[s].size();
    size_t rank = allcomm_.rank();  
    
    size_t barrier = Params::Inst()->limits.stage.barrier;

    // first nodes always gets maximum number of atoms in div assignment:
    DivAssignment zero_node_assignment(NNPP,0,NA);
    size_t maxatoms = zero_node_assignment.size();
    
    coor_t* p_alignedframe = (coor_t*) malloc((maxatoms*NNPP)*sizeof(coor_t));
    coor_t* p_alignedframeOUT = (coor_t*) malloc((maxatoms*NNPP)*sizeof(coor_t));
                    
    for(size_t f = 0; f < LNF; ++f)
    {
        // first frame = frameoffset
        size_t firstframe = framesbuffer[0][f];
        
        coor_t* p_from = &(p_coordinates_buffer[f*NA*3]);
        copyalign_frame(p_alignedframe,p_from,maxatoms);
        
        boost::mpi::all_to_all(partitioncomm_,p_alignedframe,maxatoms,p_alignedframeOUT);
        fill_alignedframe(p_alignedframeOUT,maxatoms,firstframe);
    }      
}

//void DataStagerByAtom::distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s) {
//    
//    size_t LNF = framesbuffer[s].size();
//    size_t rank = allcomm_.rank();  
//    
//    size_t barrier = Params::Inst()->limits.stage.barrier;
//                    
//    for (size_t i=0;i<NNPP;i++) {
//
//        size_t target_node = i;
//
//        if (rank==s) {
//            
//            DivAssignment target_node_assignment(NNPP,i,NA);
//            size_t off = target_node_assignment.offset();
//            size_t len = target_node_assignment.size();
//            
//            coor_t* p_fsdata = (coor_t*) malloc(LNF*len*3*sizeof(coor_t));
//            for(size_t f = 0; f < LNF; ++f)
//            {
//                coor_t* p_from = &(p_coordinates_buffer[f*NA*3+off*3]);
//                coor_t* p_to = &(p_fsdata[f*len*3]);
//                memcpy(p_to,p_from,len*3*sizeof(coor_t));
//            }
//            
//            if (target_node==s) {       
//                fill_coordinates(p_fsdata,len,framesbuffer[s]);
//            } else {
//                allcomm_.send(target_node,0,p_fsdata,LNF*len*3);
//            }
//            free(p_fsdata);
//            
//        } else if (target_node==rank) {
//            size_t len = FC_assignment.size();
//        
//            coor_t* p_localdata = (coor_t*) malloc(LNF*len*3*sizeof(coor_t));
//            allcomm_.recv(s,0,p_localdata,LNF*len*3)
//            fill_coordinates(p_localdata,len,framesbuffer[s]);
//            free(p_localdata);
//        }
//        
//        if (((i+1)%barrier)==0) {
//            timer_.start("st:dist:wait1");
//            allcomm_.barrier();
//            timer_.stop("st:dist:wait1");        
//        }
//    }
//
//    timer_.start("st:dist:wait2");
//    allcomm_.barrier();
//    timer_.stop("st:dist:wait2");        
//}

//void DataStagerByAtom::fill_coordinates(coor_t* p_localdata,size_t len,vector<size_t> frames) {
//    for(size_t f = 0; f < frames.size(); ++f)
//    {
//        for(size_t n = 0; n < len; ++n)
//        {
//            p_coordinates[ NF*n*3 + frames[f]*3 ] = p_localdata[ len*f*3 + n*3 ];
//            p_coordinates[ NF*n*3 + frames[f]*3 + 1 ] = p_localdata[ len*f*3 + n*3 + 1 ];
//            p_coordinates[ NF*n*3 + frames[f]*3 + 2 ] = p_localdata[ len*f*3 + n*3 + 2 ];
//        }
//    }
//}

void DataStagerByAtom::fill_coordinates(coor_t* p_localdata,size_t len,vector<size_t> frames) {
    for(size_t f = 0; f < frames.size(); ++f)
    {
        for(size_t n = 0; n < len; ++n)
        {
            //p_coordinates[ NF*n*3 + frames[f]*3 ] = p_localdata[ len*f*3 + n*3 ];
            //p_coordinates[ NF*n*3 + frames[f]*3 + 1 ] = p_localdata[ len*f*3 + n*3 + 1 ];
            //p_coordinates[ NF*n*3 + frames[f]*3 + 2 ] = p_localdata[ len*f*3 + n*3 + 2 ];
        }
    }
}


void DataStagerByAtom::stage_fillpartitions() {
    
    bool firstpartition=true;
    if (allcomm_.rank()>=partitioncomm_.size()) firstpartition=false;
    size_t partitions = allcomm_.size()/partitioncomm_.size();

    if (firstpartition) {
        for(size_t i = 1; i < partitions; ++i)
        {
            size_t target_node = i*partitioncomm_.size() + partitioncomm_.rank();
            allcomm_.send(target_node,0,p_coordinates,FC_assignment.size()*NF*3);            
        }
    } else {
        size_t source_node = partitioncomm_.rank();
        allcomm_.recv(source_node,0,p_coordinates,FC_assignment.size()*NF*3);            
    }
}

// end of file
