/*
 *  scatter_spectrum_writer.hpp
 *
 *  Created on: May 26, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SCATTERSPECTRUMWRITER_HPP_
#define SCATTERSPECTRUMWRITER_HPP_

// common header
#include "common.hpp"

// standard header
#include <complex>
#include <map>
#include <string>
#include <vector>

// special library headers

// other headers
#include "coor3d.hpp"

class ScatterSpectrumWriter {
	p_scatspec;
	std::string filename;
	std::string format;
public:
	ScatterSpectrumWriter(ScatterSpectrum& scatspec, std::string filename, std::string format);	

	virtual void write() = 0;
};

class ScatterSpectrumWriterPlain : ScatterSpectrumWriter {

public: 
	
	void write();
};

class ScatterSpectrumWriterFrameAverage : ScatterSpectrumWriter {

public: 
	
	void write();
};


#endif

// end of file
