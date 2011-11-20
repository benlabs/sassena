/** \file
This file contains an Interface class for communicating data structures across the network.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef MPI__WRAPPER_HPP_
#define MPI__WRAPPER_HPP_

#include "boost/mpi.hpp"

namespace mpi {
    namespace wrapper {
		/** 
		Communicates a string object across network via a MPI broadcast.
		*/
        void broadcast_stream(boost::mpi::communicator& comm,std::stringstream& stream, size_t root);        

		/** 
		Communicates a class object across network via a MPI broadcast. Class must be serializable.
		*/        
        template <class T> void broadcast_class(boost::mpi::communicator& comm,T& any, size_t root);
    }
}

#endif

// end of file