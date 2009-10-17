/*
 *  parameters.hpp
 *reg
 *  Created on: June 25, 2009
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

// other headers
#include "coor3d.hpp"

// This is a wrapper class to interface the settings implementation. The rational is to 
// move all possible configuration errors towards the initialization of the software
// preferably the 'parameters' class checks for all required settings and implementents
// default values
// Also, use hardwired constant names to move possible errors to compile time.
// Basically this class maps the structure of the configuration file, more or less

// these constructs are to be used w/ in the code the following way:
// string fs = Params::Inst()->sample.structure.file

// resolve dependency:
class CartesianCoor3D;


class SampleStructureParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & file;
		ar & format;
    }
	/////////////////// 

public:	
	std::string file;
	std::string format;
};

class SampleGroupParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & name;
		ar & file;
		ar & format;
		ar & select;
		ar & select_value;
    }
	/////////////////// 

public:	
	std::string name;
	std::string file;
	std::string format;
	std::string select;
	double select_value;
	
	SampleGroupParameters() {}
	SampleGroupParameters(const std::string& v1,const std::string& v2,const std::string& v3,const std::string& v4,const double& v5) {
		name = v1; file = v2; format = v3; select = v4; select_value = v5; 
	}
};

class SampleFramesetParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & first;
		ar & last;
		ar & last_set;
		ar & stride;
		ar & filename;
		ar & type;
    }
	/////////////////// 

public:	
	size_t first;
	size_t last;
	bool last_set;
	size_t stride;
	std::string filename;
	std::string type;
};


class SampleFramesParameters : public std::vector<SampleFramesetParameters > {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<SampleFramesetParameters>, SampleFramesParameters>(*this);
    }
	/////////////////// 

public:	
};

class SamplePBCParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & wrapping;
		ar & selection;
    }
	/////////////////// 

public:	
	bool wrapping;
	std::string selection;
};

class SampleAlignmentOriginParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & selection;
        ar & basevector;
    }
	/////////////////// 

public:	
    std::string type;
    std::string selection;
    CartesianCoor3D basevector;
};

class SampleAlignmentAxisParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & selection;
        ar & basevector;
    }
	/////////////////// 

public:	
    std::string type;
    std::string selection;
    CartesianCoor3D basevector;
};
class SampleAlignmentParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & origin;
		ar & axis;
    }
	/////////////////// 

public:	
    SampleAlignmentOriginParameters origin;
    SampleAlignmentAxisParameters axis;
};


class SampleMotionParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & displace;
		ar & direction;
		ar & frequency;
		ar & selection;
		ar & seed;
    }
	/////////////////// 

public:	
	std::string type;
	double displace;
	double frequency;
	std::string selection;
	long seed;
	CartesianCoor3D direction;
};

class SampleParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & structure;
		ar & groups;
		ar & frames;
		ar & pbc;
		ar & deuter;
		ar & motions;
        ar & alignment;
    }
	/////////////////// 

public:	
	SampleStructureParameters structure;
	std::map<std::string,SampleGroupParameters> groups;
	SampleFramesParameters frames;
	SampleAlignmentParameters alignment;	
	SamplePBCParameters pbc;
	std::vector<SampleMotionParameters> motions;
	std::vector<std::string> deuter;
};

class ScatteringBackgroundPhaseParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & selection;
		ar & scaling;
		ar & factor;	
		ar & nullrange;	
    }
	/////////////////// 

public:	
	ScatteringBackgroundPhaseParameters()  {
		scaling = "auto";
		factor = 1.0;
	}

	std::string selection;
	double factor;
	std::string scaling;
	double nullrange;
};

class ScatteringBackgroundParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & factor;
		ar & phases;
		ar & gridspacing;
		ar & stride;
    }
	/////////////////// 

public:	
	ScatteringBackgroundParameters() {
		stride = 1;
		gridspacing = 3;
		factor = 0.0;
		type = "manual";
	}

	std::string type;
	double factor;

	std::vector<ScatteringBackgroundPhaseParameters> phases;
	
	double gridspacing;
	size_t stride;
};


class ScatteringAverageOrientationParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & method;
		ar & vectors;
		ar & resolution;
		ar & origin;
		ar & axis;
    }
	/////////////////// 

public:	
	ScatteringAverageOrientationParameters() {
		type = "none";
		method = "bruteforce";
		vectors = "mcboostunisphere";
		resolution = 1000;
		origin = "";
		axis = CartesianCoor3D(1,0,0);
	}
	
	std::string type;
	std::string method;
	std::string vectors;
	double resolution;
	
	std::string origin;
	CartesianCoor3D axis;
	
};

class ScatteringAverageParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & orientation;
    }
	/////////////////// 

public:	
	ScatteringAverageOrientationParameters orientation;
};

class ScatteringCorrelationParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & method;
        ar & zeromean;
    }
	/////////////////// 

public:	
	ScatteringCorrelationParameters() {
		type = "none";
		method = "direct";
	}
	std::string type;
	std::string method;
    bool zeromean;	
};

class ScatteringInterferenceParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
    }
	/////////////////// 

public:	
	ScatteringInterferenceParameters() {
		type = "all";
	}
	std::string type;
};

class ScatteringVectorsScanParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & from;
		ar & to;
		ar & points;
		ar & exponent;
		ar & basevector;	
    }
	/////////////////// 

public:	
	ScatteringVectorsScanParameters() : points(100), exponent(1.0), basevector(CartesianCoor3D(1,0,0)) {}
	
	double from;
	double to;
	size_t points;
	double exponent;
	CartesianCoor3D basevector;
};

class ScatteringVectorsParameters : public std::vector<CartesianCoor3D> {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<CartesianCoor3D> >(*this);
		ar & scans;
    }
	/////////////////// 

public:
	std::vector<ScatteringVectorsScanParameters> scans;
	
	void create_from_scans();
};

class ScatteringParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & scatterfactors;
		ar & interference;
		ar & correlation;
		ar & target;
		ar & qvectors;
		ar & average;
		ar & background;
    }
	/////////////////// 

public:
	
	ScatteringInterferenceParameters interference;
	ScatteringCorrelationParameters correlation;
	std::string target;

	std::string scatterfactors;
	
	ScatteringVectorsParameters qvectors;
	
	ScatteringAverageParameters average;
	ScatteringBackgroundParameters background;	
};

class OutputFileParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & method;
		ar & format;
		ar & filename;
		ar & name;
    }
	/////////////////// 

public:
	
	std::string name;
	std::string method;
	std::string format;
	std::string filename;
};

class OutputParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & prefix;
		ar & files;
    }
	/////////////////// 

public:
	OutputParameters() {
		prefix = "scattering-data-";
	}
	
	std::string prefix;
	std::vector<OutputFileParameters> files;
};

class LimitsBuffersParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & allgather_max;
    }
	/////////////////// 

public:
	size_t allgather_max;
};

class LimitsParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & framecache_max;
		ar & coordinatesets_cache_max;		
		ar & static_load_imbalance_max;
		ar & buffers;
    }
	/////////////////// 

public:
	size_t framecache_max;
	size_t coordinatesets_cache_max;
	double static_load_imbalance_max;
	LimitsBuffersParameters buffers;
};


// this section is dedicated to parameters which are computed on the fly!
class RuntimeParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & config_rootpath;
    }
	/////////////////// 

public:
	std::string config_rootpath;
};

class DebugParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & timer;
		ar & barriers;
		ar & scatter_from_frame;
		ar & decomposition_split;
    }
	/////////////////// 

public:
	bool timer;
	bool barriers;
	bool scatter_from_frame;
	size_t decomposition_split;
};




// implement Parameters as a singleton class, this makes it globally available
// it requires the call of the init() function, which also implements all checks
class Params  {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & carboncopy;
		ar & config_rootpath;
		ar & sample;
		ar & runtime;
		ar & scattering;
		ar & output;
		ar & limits;
		ar & debug;
    }
	/////////////////// 

	Params() {}
	Params(const Params&);
	Params& operator=(const Params&);

	std::vector<std::string> carboncopy;
	
	std::string config_rootpath;
	
	std::string get_filepath(std::string filename);	
	
	void read_conf(std::string filename);
	void read_xml(std::string filename);

	std::string guessformat(std::string filename);
	bool check();
	
public: 
	// interface for parameters
	RuntimeParameters runtime;
	SampleParameters sample;
	ScatteringParameters scattering;
	OutputParameters output;
	LimitsParameters limits;
	DebugParameters debug;
	
	// interface for initiatilzation and interfacing
	static Params* Inst() { static Params instance; return &instance;}
	
	void init(std::string filename);
	~Params() {}; // it is said some compilers have problems w/ private destructors.
	
	void write(std::string filename, std::string format="");	
};

#endif 

// end of file
