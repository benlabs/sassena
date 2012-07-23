/** \file
This file contains a class which implements the data staging logic for the trajectory data. The first partition reads the trajectory data from disc and distributes the trajectory to any other partition. 

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "stager/data_stager.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/thread.hpp>

// other headers
#include "stager/coordinate_writer.hpp"
#include "math/coor3d.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;

/** 
Initializes a data staging object which decomposes the trajectory by frames using div logic (consecutive frames are assigned to the same node)
\param[in] sample Sample interface, which is used to load the frames.
\param[in] allcomm MPI Communicator which contains all nodes which participate in the analysis
\param[in] partitioncomm MPI Communicator which contains all nodes of the local partition.
\param[in,out] timer Timer object which is used to store timing information 
*/
DataStagerByFrame::DataStagerByFrame(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,Timer& timer) 
    : m_sample(sample),
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),    
    timer_(timer)
{
    NN = allcomm_.size();  
    NNPP = partitioncomm_.size();
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->stager.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = allcomm_.rank();

	DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);

    // check data size requirements & allocate memory
	size_t data_bytesize = assignment.max()*NA*3*sizeof(coor_t);
    if (Params::Inst()->limits.stage.memory.data<data_bytesize) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize));
        }
        throw;
    }

    p_coordinates = (coor_t*) malloc(data_bytesize);
}

/** 
Triggers the staging procedure. Data is staged on the first partition and then cloned to all other partitions. 
\return Pointer to the local memory which contains the staged coordinates (2D data, atoms x assigned frames)
*/
coor_t* DataStagerByFrame::stage() {
    
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

    if (Params::Inst()->stager.dump) {
        Info::Inst()->write(string("Dumping coordinates to file: ")+Params::Inst()->stager.filepath);
        write(Params::Inst()->stager.filepath,Params::Inst()->stager.format);
    }

    return p_coordinates;
}

/** 
Loads trajectory data into the first partition.
*/
void DataStagerByFrame::stage_firstpartition() {
    
	DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);
    
    for(size_t f=0;f<assignment.size();f++) {
        timer_.start("st:load");
        CoordinateSet* p_cset = m_sample.coordinate_sets.load(assignment[f]);
                    
        for(size_t n=0;n<NA;n++) {
            p_coordinates[f*3*NA+3*n]=p_cset->c1[n];
            p_coordinates[f*3*NA+3*n+1]=p_cset->c2[n];
            p_coordinates[f*3*NA+3*n+2]=p_cset->c3[n];            
        }
        timer_.stop("st:load");
        delete p_cset;
    }
}

/** 
Clones coordinates from the first to all remaining partitions.
*/
void DataStagerByFrame::stage_fillpartitions() {
	DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);

    // create a communicator to broadcast between partitions.    
    boost::mpi::communicator interpartitioncomm_ = allcomm_.split(partitioncomm_.rank());
    boost::mpi::broadcast(interpartitioncomm_,p_coordinates,assignment.size()*NA*3,0);
}

/** 
Dumps the coordinates to an output file. Writes in parallel.
*/
void DataStagerByFrame::write(std::string filename, std::string format) {
	DivAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NF);

    // use the first partition by default
    if (allcomm_.rank()<partitioncomm_.size()) {
        timer_.start("st:dump");

        ICoordinateWriter* p_cw = NULL;
        if (format=="dcd") {
            p_cw = new DCDCoordinateWriter(filename,NF,NA);
        } else {
            Err::Inst()->write(string("Format for coordinate dumping not known: ")+format);
            throw;
        }
        if (allcomm_.rank()==0) p_cw->init();
        partitioncomm_.barrier();
        p_cw->prepare();
        partitioncomm_.barrier();

		// synchronize writes, to write a consecutive series of blocks in parallel
		for(size_t c = 0; c < assignment.max(); ++c)
		{
		   if (c<assignment.size()) {
				p_cw->write(&(p_coordinates[c*NA*3]),assignment[c],1);		
		   }
		   partitioncomm_.barrier();
		}
        timer_.stop("st:dump");
    }

    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");
}

/** 
Initializes a data staging object which decomposes the trajectory by atoms using modulo logic (consecutive atoms are assigned to different nodes)
\param[in] sample Sample interface, which is used to load the frames.
\param[in] allcomm MPI Communicator which contains all nodes which participate in the analysis
\param[in] partitioncomm MPI Communicator which contains all nodes of the local partition.
\param[in,out] timer Timer object which is used to store timing information 
*/
DataStagerByAtom::DataStagerByAtom(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,Timer& timer) 
    : m_sample(sample),
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),    
    timer_(timer)
{
    NN = allcomm_.size(); 
    NNPP = partitioncomm_.size();    
    NP = NN/NNPP;
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->stager.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = allcomm_.rank();

	// assignment for atom decomposition is modulo
	ModAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NA);

    // check data size requirements & allocate memory
    size_t data_bytesize = assignment.size()*NF*3*sizeof(coor_t);
    size_t data_bytesize_indicator_max = 0;
    boost::mpi::all_reduce(partitioncomm_,data_bytesize,data_bytesize_indicator_max,boost::mpi::maximum<size_t>());

    if (Params::Inst()->limits.stage.memory.data<data_bytesize_indicator_max) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
        }
        throw;
    }

    p_coordinates = (coor_t*) malloc(data_bytesize);
}

/** 
Triggers the staging procedure. Data is staged on the first partition and then cloned to all other partitions. 
\return Pointer to the local memory which contains the staged coordinates (2D data, frames x assigned atoms)
*/
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

    if (Params::Inst()->stager.dump) {
	    if (allcomm_.rank()==0) {
	        Info::Inst()->write(string("Dumping coordinates to file: ")+Params::Inst()->stager.filepath);		
		}
        write(Params::Inst()->stager.filepath,Params::Inst()->stager.format);
    }

	ModAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NA);


    return p_coordinates;
}

/** 
Loads trajectory data into the first partition. Performs transposition of the data on the fly by making heavy use of buffers and MPI all to all communication.
*/
void DataStagerByAtom::stage_firstpartition() {
    // determine frame assignment (for reading) for this node
	DivAssignment frameassignment(partitioncomm_.size(),partitioncomm_.rank(),NF);

	// test buffer size, minimum allowed is 1 frame
    size_t buffer_bytesize = Params::Inst()->limits.stage.memory.buffer;
    size_t frame_bytesize = NA*3*sizeof(coor_t);
	size_t framebuffer_max = buffer_bytesize/frame_bytesize;

    if (framebuffer_max==0) {
        if (partitioncomm_.rank()==0) {
            Err::Inst()->write("Cannot load trajectory into buffer.");
            Err::Inst()->write(string("limits.memory.data_stager=")+boost::lexical_cast<string>(Params::Inst()->limits.stage.memory.buffer));
            Err::Inst()->write(string("requested=")+boost::lexical_cast<string>(frame_bytesize));            
        }
        throw;
    }
		
	// determine whether maxframes_read_total fits into the buffer
	// if not: we need to perform multiple cycles
	size_t cycles=1;
 	if (frameassignment.max()>framebuffer_max) {
		cycles=frameassignment.max()/framebuffer_max;
		if ((frameassignment.max()%framebuffer_max)!=0) cycles+=1;
	} else {
		// reduce the buffer size to the maximum required
		framebuffer_max = frameassignment.max();
	}

	// allocate buffer
    if (partitioncomm_.rank()==0) {
        Info::Inst()->write(string("Initializing buffer size to: ")+boost::lexical_cast<string>(framebuffer_max));
    }
	// align the readbuffer with the partition size, makes sure that access to the memory is valid
	size_t NA_aligned = NA;
	if ( (NA_aligned%NNPP)!=0 ) {
		NA_aligned = ((NA/NNPP) +1 )*NNPP;
	}
    coor_t* p_coordinates_readbuffer = (coor_t*) malloc(framebuffer_max*NA_aligned*3*sizeof(coor_t));
    coor_t* p_coordinates_exchangebuffer = (coor_t*) malloc(framebuffer_max*NNPP*3*sizeof(coor_t));

	for (size_t c=0;c<cycles;c++) {
		// fill buffer
		for (size_t f=framebuffer_max*c;f<std::min(framebuffer_max*(c+1),frameassignment.size());f++) {
			timer_.start("st:load");
		    CoordinateSet* p_cset = m_sample.coordinate_sets.load(frameassignment[f]);
			size_t framepos_buffer = f-framebuffer_max*c;
		    for(size_t n=0;n<NA;n++) {
				// write it transposed, makes exchange through MPI easy!
		        p_coordinates_readbuffer[n*framebuffer_max*3 + framepos_buffer*3   ]=p_cset->c1[n];
		        p_coordinates_readbuffer[n*framebuffer_max*3 + framepos_buffer*3 +1]=p_cset->c2[n];
		        p_coordinates_readbuffer[n*framebuffer_max*3 + framepos_buffer*3 +2]=p_cset->c3[n];            
		    }
		    delete p_cset;
		    timer_.stop("st:load");		
		}
		// exchange data with other nodes
		for (size_t ec=0;ec<(NA_aligned/NNPP);ec++) {
			// places atoms in modulo fashion on different nodes
			boost::mpi::all_to_all(
				partitioncomm_, 
				&( p_coordinates_readbuffer[ ec*framebuffer_max*3*NNPP ] ), 
				framebuffer_max*3, 
				p_coordinates_exchangebuffer
			);			
			// data is now in exchangebuffer,
			// organized as time sequence fragments for a particular atom (id=NNPP*ec+partitioncomm_.rank())
			size_t atomid=NNPP*ec+partitioncomm_.rank();
			if (atomid<NA) {
				//iterate over all time fragments
				for (size_t n=0;n<NNPP;n++) {
					// get properties of the time fragment from associated assignment
					DivAssignment nodeframeassignment(NNPP,n,NF);
					size_t len = nodeframeassignment.size();
					if ( framebuffer_max<nodeframeassignment.size() ) {
						len = std::min(framebuffer_max,nodeframeassignment.size()-c*framebuffer_max);
					}
					size_t pos = nodeframeassignment.offset()+c*framebuffer_max;
					
					coor_t* p_to = &( p_coordinates[ (ec*NF + pos )*3 ] );
					coor_t* p_from = &( p_coordinates_exchangebuffer[ n*framebuffer_max*3 ] );
					memcpy(p_to,p_from,3*len*sizeof(coor_t));					
				}
			}
		}
	}
	
	free(p_coordinates_readbuffer);
	free(p_coordinates_exchangebuffer);
}

/** 
Clones coordinates from the first to all remaining partitions.
*/
void DataStagerByAtom::stage_fillpartitions() {   
	ModAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NA);
	
    // create a communicator to broadcast between partitions.    
    boost::mpi::communicator interpartitioncomm_ = allcomm_.split(partitioncomm_.rank());
    boost::mpi::broadcast(interpartitioncomm_,p_coordinates,assignment.size()*NF*3,0);
}

/** 
Dumps the coordinates to an output file. Writes in parallel.
*/
void DataStagerByAtom::write(std::string filename, std::string format) {
	ModAssignment assignment(partitioncomm_.size(),partitioncomm_.rank(),NA);
	
    // use the first partition by default
    if (allcomm_.rank()<partitioncomm_.size()) {
        timer_.start("st:dump");

        ICoordinateWriter* p_cw = NULL;
        if (format=="dcd") {
            p_cw = new DCDCoordinateWriter(filename,NA,NF);
        } else {
			if (allcomm_.rank()==0) {
	            Err::Inst()->write(string("Format for coordinate dumping not known: ")+format);				
			}
            throw;
        }
        if (partitioncomm_.rank()==0) p_cw->init();
        partitioncomm_.barrier();
        p_cw->prepare();
        partitioncomm_.barrier();

		// synchronize writes, to write a consecutive series of blocks in parallel
		for(size_t c = 0; c < assignment.max(); ++c)
		{
		   if (c<assignment.size()) {
				p_cw->write(&(p_coordinates[c*NF*3]),assignment[c],1);		
		   }
		   partitioncomm_.barrier();
		}
		timer_.stop("st:dump");
    }

    timer_.start("st:wait");
    allcomm_.barrier();
    timer_.stop("st:wait");
}

// end of file
