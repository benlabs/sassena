/** \file
This file contains the parameters class which contains the settings and values used to adjust the program control flow to the user input files.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef CONTROL__PARAMETERS_HPP_
#define CONTROL__PARAMETERS_HPP_

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
#include <boost/program_options.hpp>

// other headers
#include "math/coor3d.hpp"

// This class follow a strong hierarchy. The root class is located at the end of this file.

// resolve dependency:
class CartesianCoor3D;

/**
Section which defines the structure
*/
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

/**
Section which defines a generic selection (used as parent by specific selections)
*/
class SampleSelectionParameters {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type_;
    }
    std::string type_;
public:
    SampleSelectionParameters() {}
    SampleSelectionParameters(std::string type) : type_(type) {}
    
    std::string type() { return type_;}
};

/**
Section which defines a selection based on individual indexes
*/
class SampleIndexSelectionParameters : public SampleSelectionParameters {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<SampleSelectionParameters, SampleIndexSelectionParameters>(*this);
		ar & ids_;
    }
public:
    
    std::vector<size_t> ids_;
    SampleIndexSelectionParameters() {}    
    SampleIndexSelectionParameters(std::vector<size_t> ids) : SampleSelectionParameters("index"), ids_(ids) {}
};

/**
Section which defines a selection based on a given range
*/
class SampleRangeSelectionParameters : public SampleSelectionParameters {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<SampleSelectionParameters, SampleRangeSelectionParameters>(*this);
		ar & from_;
        ar & to_;
    }
public:
    
    size_t from_;
    size_t to_;

    SampleRangeSelectionParameters() {}
    SampleRangeSelectionParameters(size_t from, size_t to) : SampleSelectionParameters("range"), from_(from), to_(to) {}
};

/**
Section which defines a selection based on a lexical pattern (regular expression matching atom labels)
*/
class SampleLexicalSelectionParameters : public SampleSelectionParameters {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<SampleSelectionParameters, SampleLexicalSelectionParameters>(*this);
		ar & expression_;
    }
public:
    
    std::string expression_;
    
    SampleLexicalSelectionParameters() {}
    SampleLexicalSelectionParameters(std::string expression) : SampleSelectionParameters("lexical"), expression_(expression) {}
};

/**
Section which defines a file based selection
*/
class SampleFileSelectionParameters : public SampleSelectionParameters {
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<SampleSelectionParameters, SampleFileSelectionParameters>(*this);
		ar & file_;
        ar & format_;
        ar & selector_;
        ar & expression_;
    }
public:

    std::string file_;
    std::string format_;
    std::string selector_;    
    std::string expression_;
    
    SampleFileSelectionParameters() {}
    SampleFileSelectionParameters(std::string file, std::string format, std::string selector, std::string expression) : SampleSelectionParameters("file"), file_(file), format_(format), selector_(selector), expression_(expression) {}
};

/**
Section which defines a single trajectory file entry
*/
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
        ar & clones;
        ar & index;
		ar & index_default;
    }
	/////////////////// 

public:	
	size_t first;
	size_t last;
    size_t clones;
	bool last_set;
	size_t stride;
	std::string filename;
	std::string type;
    std::string index;
	bool index_default;
};

/**
Section which lists the used trajectory files
*/
class SampleFramesetsParameters : public std::vector<SampleFramesetParameters > {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<SampleFramesetParameters>, SampleFramesetsParameters>(*this);
    }
	/////////////////// 

public:	
};

/**
Section which stores reference information which may be required during some alignment procedures
*/
class SampleMotionReferenceParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & frame;
        ar & file;
        ar & format;
        ar & selection;
    }
	/////////////////// 

public:	
	std::string type;
	std::string selection;
	std::string file;
	std::string format;
	size_t frame;
};
/**
Section which defines artificial motions
*/
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
        ar & radius;
		ar & selection;
		ar & seed;
        ar & sampling;
		ar & reference;
    }
	/////////////////// 

public:	
	std::string type;
	double displace;
	double frequency;
    double radius;
	std::string selection;
	long seed;
	long sampling;	
	CartesianCoor3D direction;
    SampleMotionReferenceParameters reference;
};

/**
Section which stores reference information which may be required during some alignment procedures
*/
class SampleAlignmentReferenceParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & frame;
        ar & file;
        ar & format;
        ar & selection;
    }
	/////////////////// 

public:	
	std::string type;
	std::string selection;
	std::string file;
	std::string format;
	size_t frame;
};

/**
Section which stores alignment information, which is applied during the staging of the trajectory data
*/
class SampleAlignmentParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & selection;
        ar & order;
        ar & reference;
    }
	/////////////////// 

public:	
    SampleAlignmentReferenceParameters reference;
	std::string type;
	std::string selection;
    std::string order;
};

/**
Section which stores sample specific parameters
*/
class SampleParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<SampleIndexSelectionParameters*>(NULL));
        ar.register_type(static_cast<SampleRangeSelectionParameters*>(NULL));
        ar.register_type(static_cast<SampleFileSelectionParameters*>(NULL));
        ar.register_type(static_cast<SampleLexicalSelectionParameters*>(NULL));
        
		ar & structure;
		ar & selections;
		ar & framesets;
		ar & motions;
        ar & alignments;
    }
	/////////////////// 

public:	
    ~SampleParameters();
    
	SampleStructureParameters structure;
	std::map<std::string,SampleSelectionParameters*> selections;
	SampleFramesetsParameters framesets;
	std::vector<SampleMotionParameters> motions;
    std::vector<SampleAlignmentParameters> alignments;
};

/**
Section which stores selection based scaling factors for background correction
*/
class ScatteringBackgroundKappaParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & selection;
        ar & value;
    }
	/////////////////// 

public:	
	std::string selection;
	double value;
};

/**
Section which stores background correction parameters
*/
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
		ar & kappas;
    }
	/////////////////// 

public:	
    std::string type;
    double factor;
    
    std::vector<ScatteringBackgroundKappaParameters> kappas;
};

/**
Section which is used when vector based orientational averaging is performed
*/
class ScatteringAverageOrientationVectorsParameters : public std::vector<CartesianCoor3D> {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object<std::vector<CartesianCoor3D> >(*this);
		ar & type;
		ar & algorithm;
		ar & resolution;
		ar & file;
		ar & seed;
    }
	/////////////////// 

public:	
	std::string type;
	std::string algorithm;
	std::string file;
	size_t resolution;
	long seed;
	
	void create();
};

/**
Section which is used to store the used multipole identifiers when multipole based orientational averaging is performed
*/
class ScatteringAverageOrientationMultipoleMomentsParameters : public std::vector<std::pair<long,long> > {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & resolution;
        ar & file;
    }
	/////////////////// 

public:	
	std::string type;
    long resolution;
    std::string file;

	void create();
};

/**
Section which is used when multipole based orientational averaging is performed
*/
class ScatteringAverageOrientationMultipoleParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & moments;
    }
	/////////////////// 

public:	
	std::string type;
	ScatteringAverageOrientationMultipoleMomentsParameters moments;
};

/**
Section which is used when exact orientational averaging is performed
*/
class ScatteringAverageOrientationExactParameters {
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
	std::string type; 
};

/**
Section which defines orientational averaging procedures
*/
class ScatteringAverageOrientationParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
        ar & axis;

		ar & vectors;
		ar & multipole;
		ar & exact;
    }
	/////////////////// 
	
public:		
	std::string type; 
	CartesianCoor3D axis;
	
	ScatteringAverageOrientationVectorsParameters vectors;
	ScatteringAverageOrientationMultipoleParameters multipole;
	ScatteringAverageOrientationExactParameters exact;

};

/**
Section which defines averaging procedures
*/
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

/**
Section which defines further processing of the scattering signal, e.g. autocorrelation
*/
class ScatteringDspParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & method;
    }
	/////////////////// 

public:	
	std::string type;
	std::string method;
};

/**
Section which describes a scan entry used to construct the scattering vectors
*/
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

/**
Section which defines the scattering vectors
*/
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

/**
Section which stores parameters used during the writing of the signal file
*/
class ScatteringSignalParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & file;
        ar & fqt;
        ar & fq;
        ar & fq0;
        ar & fq2;
    }
	/////////////////// 

public:
	std::string file;
    bool fqt;
    bool fq;
    bool fq0;
    bool fq2;
};

/**
Section which stores parameters used during the scattering calculation
*/
class ScatteringParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & type;
		ar & dsp;
		ar & qvectors;
		ar & average;
		ar & background;
        ar & signal;
    }
	/////////////////// 

public:	
	std::string type;
	ScatteringDspParameters dsp;

	ScatteringVectorsParameters qvectors;
	
	ScatteringAverageParameters average;
	ScatteringBackgroundParameters background;	

	ScatteringSignalParameters signal;	
};

/**
Section which stores parameters used during data staging
*/
class StagerParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & dump;
		ar & file;
		ar & format;
		ar & target;
        ar & mode;
    }
	/////////////////// 

public:	
    bool dump;
	std::string file;
	std::string format;
	std::string target;
    std::string mode;
};

/**
Section which stores parameters affecting the memory limitations of the file writer service
*/
class LimitsServicesSignalMemoryParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & server;
        ar & client;
    }
	/////////////////// 

public:
    size_t server;
    size_t client;
};

/**
Section which stores parameters affecting the timing of the file writer service
*/
class LimitsServicesSignalTimesParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & serverflush;
        ar & clientflush;
    }
	/////////////////// 

public:
    size_t serverflush;
    size_t clientflush;
};

/**
Section which stores parameters affecting the file writer service, which writes results to the signal file
*/
class LimitsServicesSignalParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & memory;
        ar & times;
    }
	/////////////////// 

public:
    LimitsServicesSignalMemoryParameters memory;
    LimitsServicesSignalTimesParameters times;    
};

/**
Section which stores parameters affecting the monitoring service, which reports progress to the console
*/
class LimitsServicesMonitorParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & delay;
        ar & sampling;
    }
	/////////////////// 

public:
    size_t delay;
    size_t sampling;
};

/**
Section which stores parameters affecting the services
*/
class LimitsServicesParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & signal;
        ar & monitor;
    }
	/////////////////// 

public:
    LimitsServicesSignalParameters signal;
    LimitsServicesMonitorParameters monitor;
};

/**
Section which stores memory limits during the computation
*/
class LimitsComputationMemoryParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & result_buffer;
        ar & alignpad_buffer;
        ar & exchange_buffer;
		ar & signal_buffer;
		ar & scale;
    }
	/////////////////// 

public:
    size_t result_buffer;
    size_t alignpad_buffer;
	size_t exchange_buffer;
	size_t signal_buffer;
	size_t scale;
};

/**
Section which stores parameters used during the computation
*/
class LimitsComputationParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & threads;
		ar & memory;
        ar & processes;
        ar & cores;
    }
	/////////////////// 

public:
    size_t threads;
    size_t processes;
    size_t cores;
    LimitsComputationMemoryParameters memory;
};

/**
Section which stores parameters determining the computational partition size
*/
class LimitsDecompositionPartitionsParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & automatic;        
        ar & size;
    }
	/////////////////// 

public:
    bool automatic;
    size_t size;
    
};

/**
Section which stores parameters used for partitioning the computation among the available compute nodes
*/
class LimitsDecompositionParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & utilization;
        ar & partitions;        
    }
	/////////////////// 

public:
    double utilization;
    LimitsDecompositionPartitionsParameters partitions;
};

/**
Section which stores parameters used during the writing of the signal output file
*/
class LimitsSignalParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & chunksize;  
    }
	/////////////////// 

public:
    size_t chunksize;
};

/**
Section which stores memory limits during the data staging process.
*/
class LimitsStageMemoryParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & data;
        ar & buffer;
    }
	/////////////////// 

public:
    size_t data;
    size_t buffer;
};

/**
Section which stores limits regarding the data staging process.
*/
class LimitsStageParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & memory;
    }
	/////////////////// 

public:
    LimitsStageMemoryParameters memory;
};

/**
Section which stores computational limitations and performance characteristics.
*/
class LimitsParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & stage;
        ar & signal;
        ar & services;
        ar & decomposition;
        ar & computation;
    }
	/////////////////// 

public:
    LimitsStageParameters stage;    
    LimitsSignalParameters signal;    
    LimitsServicesParameters services;    
    LimitsComputationParameters computation;    
    LimitsDecompositionParameters decomposition;
};

/**
Section which stores parameters influencing the progress monitoring
*/
class DebugMonitorParameters {
   private:
   	/////////////////// MPI related
   	// make this class serializable to 
   	// allow sample to be transmitted via MPI
       friend class boost::serialization::access;	
   	template<class Archive> void serialize(Archive & ar, const unsigned int version)
       {
   		ar & update;
       }
   	/////////////////// 

   public:
   	bool update; 
};

/**
Section which stores switches for dumping information to console output.
*/
class DebugPrintParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & orientations;
    }
	/////////////////// 

public:
	bool orientations;
};

/**
Section which stores IO write specific debug parameters. Used to tune the frequency by which the results are written to the signal file. Can also be used to avoid writing to the signal file. 
*/
class DebugIowriteParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & write;
        ar & buffer;
    }
	/////////////////// 

public:
	bool write;
    bool buffer;
};

/**
Section which stores debug parameters
*/
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
        ar & monitor;
        ar & iowrite;
        ar & print;
    }
	/////////////////// 

public:
	bool timer;
	bool barriers;
    DebugIowriteParameters iowrite;
    DebugPrintParameters print;
    DebugMonitorParameters monitor;
};

/**
Section which stores a reference to the used database
*/
class DatabaseParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & type;
		ar & file;
        ar & format;
    }
	/////////////////// 

public:
	std::string type;
	std::string file;
	std::string format;
};


/**
This is a wrapper class to interface the settings implementation. The rational is to 
move all possible configuration errors towards the initialization of the software
preferably the 'parameters' class checks for all required settings and implementents
default values
Also, use hardwired constant names to move possible errors to compile time.
Basically this class maps the structure of the configuration file, more or less

these constructs are to be used w/ in the code the following way:
\verbatim
string fs = Params::Inst()->sample.structure.file
*/

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
		ar & scattering;
        ar & stager;
        ar & database;
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

	void read_xml(std::string filename);
    boost::program_options::options_description options();
    void overwrite_options(boost::program_options::variables_map& vm);
public: 
	// interface for parameters
	SampleParameters sample;
	ScatteringParameters scattering;
    StagerParameters stager;
    DatabaseParameters database;
	LimitsParameters limits;
	DebugParameters debug;
	
	// interface for initiatilzation and interfacing
	static Params* Inst() { static Params instance; return &instance;}
	
	void init(int argc,char** argv);
	~Params() {}; // it is said some compilers have problems w/ private destructors.
	
	void write_xml(std::string filename);	
};

#endif 

// end of file
