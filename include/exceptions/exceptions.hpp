/** \file
This file contains exceptions defines within the context of the application.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef EXCEPTIONS__EXCEPTIONS_HPP_
#define EXCEPTIONS__EXCEPTIONS_HPP_

// common header
#include "common.hpp"

// standard header
#include <exception>

// special library headers

// other headers

namespace sassena {
	/** 
	 Signal which should be caught in the main routine to allow a clean termination.
	*/
    class terminate_request : public std::runtime_error {
    public:
        terminate_request() : std::runtime_error("Terminate requested") {}
    };
}


#endif

// end of file