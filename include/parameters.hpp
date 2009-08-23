/*
 *  scatterdata.hpp
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
#include <libconfig.h++>

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
		ar & boost::serialization::base_object<std::vector<SampleFramesetParameters> >(*this);
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
		ar & center;
    }
	/////////////////// 

public:	
	bool wrapping;
	std::string center;
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
    }
	/////////////////// 

public:	
	SampleStructureParameters structure;
	std::map<std::string,SampleGroupParameters> groups;
	SampleFramesParameters frames;
	SamplePBCParameters pbc;
	std::vector<std::string> deuter;
};

class ScatteringBackgroundParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & method;
		ar & atom_sizes;
		ar & rescale;
		ar & include;
		ar & exclude;
		ar & value;
		ar & resolution;
		ar & hydration;
		ar & framestride;
    }
	/////////////////// 

public:	
	std::string method;
	std::string atom_sizes;
	bool rescale;
	std::vector<std::string> include;
	std::vector<std::string> exclude;
	
	double value;
	
	size_t resolution;
	double hydration;
	size_t framestride;
};

class ScatteringAverageParameters {
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
	std::string type;
	std::string method;
	std::string vectors;
	double resolution;
	
	std::string origin;
	CartesianCoor3D axis;
	
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
    }
	/////////////////// 

public:	
	std::string type;
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
	std::string type;
};


class ScatteringParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & probe;
		ar & interference;
		ar & correlation;
		ar & target;
		ar & qqqvectors;
		ar & framestride;
		ar & average;
		ar & background;
    }
	/////////////////// 

public:
	ScatteringInterferenceParameters interference;
	ScatteringCorrelationParameters correlation;
	std::string target;

	std::string probe;
	
	std::vector<CartesianCoor3D> qqqvectors;
	size_t framestride;
	
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
		ar & type;
		ar & format;
		ar & filename;
    }
	/////////////////// 

public:
	std::string type;
	std::string format;
	std::string filename;
	
	OutputFileParameters() {}	
	OutputFileParameters(const std::string& v1,const std::string& v2,const std::string& v3) {
		type = v1; format = v2; filename = v3;
	}
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
		ar & static_load_imbalance_max;
		ar & buffers;
    }
	/////////////////// 

public:
	size_t framecache_max;
	double static_load_imbalance_max;
	LimitsBuffersParameters buffers;
};



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
    }
	/////////////////// 

public:
	bool timer;
	bool barriers;
	bool scatter_from_frame;
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
		ar & scattering;
		ar & output;
		ar & limits;
    }
	/////////////////// 

	Params() {}
	Params(const Params&);
	Params& operator=(const Params&);

	std::vector<std::string> carboncopy;
	
	std::string config_rootpath;
	
	std::string get_filepath(std::string filename);	
	
	void read_conf(std::string filename);

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
	
	void init(std::string filename, std::string format="");
	~Params() {}; // it is said some compilers have problems w/ private destructors.
	
	void write(std::string filename, std::string format="");	
};


class DatabaseVolumesParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & m_functiontypes;
		ar & m_constants;
		ar & m_quicklookup;
		ar & quicklookup_counter;
    }
	/////////////////// 

	std::map<size_t,size_t> m_functiontypes;
	std::map<size_t,std::vector<double> > m_constants;
	std::map<size_t,double> m_quicklookup; // use a quicklookup table for 
	size_t quicklookup_counter;

	void add_quicklookup(size_t ID, double value); 
	void clear_quicklookup();

public:	
	DatabaseVolumesParameters() : quicklookup_counter(0) {}
	
//	void reg(std::string label, std::vector<double> constants,size_t function_type);
//	double get(std::string label);
	void reg(size_t ID, std::vector<double> constants,size_t function_type);
	double get(size_t ID);
	
};

class DatabaseSFactorsParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & m_functiontypes;
		ar & m_constants;
		ar & m_quicklookup;
		ar & quicklookup_counter;
    }
	/////////////////// 
	std::map<size_t, size_t>               m_functiontypes;
	std::map<size_t, std::vector<double> > m_constants;
	
	std::map<double, std::map<size_t, double> > m_quicklookup;
	size_t quicklookup_counter;
	
	void add_quicklookup(size_t ID, double q, double value); 
	void clear_quicklookup();	
public:	
	DatabaseSFactorsParameters() : quicklookup_counter(0) {}
	
//	void reg(std::string label, std::vector<double> constants,size_t function_type);
//	double get(std::string label,double q);
	void reg(size_t ID, std::vector<double> constants,size_t function_type);
	double get(size_t ID,double q);
	
};

class DatabaseMassesParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & m_masses;
    }
	/////////////////// 

	std::map<size_t, double> m_masses;
public:	
	void reg(size_t ID, double); 	
	double get(size_t ID);		
};


class DatabaseAtomIDsParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & m_IDs;
		ar & m_rIDs;		
		ar & nextID;
    }
	/////////////////// 

	std::map<std::string,size_t> m_IDs;
	std::map<size_t,std::string> m_rIDs;
	
	size_t nextID;
public:		
	DatabaseAtomIDsParameters() : nextID(0) {}
	
	void reg(std::string label); 
	size_t get(std::string label); 	
	std::string rget(size_t ID); 
	
};

class DatabaseNamesPDBParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & m_label2regexp;
		ar & m_quicklookup;
		ar & quicklookup_counter;
    }
	/////////////////// 

	std::map<std::string,std::string> m_label2regexp;
	std::map<std::string,std::string> m_quicklookup; // use a quicklookup table for 
	
	size_t quicklookup_counter;
	
	void add_quicklookup(std::string testlabel, std::string label); 
	void clear_quicklookup();	
public:				
	DatabaseNamesPDBParameters() : quicklookup_counter(0) {}
	
	void reg(std::string label, std::string regexp); // register a regexp for a label
	std::string get(std::string testlabel); // resolve a testlabel and give label back
	
};

class DatabaseNamesParameters {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & pdb;
    }
	/////////////////// 

public:
	DatabaseNamesPDBParameters pdb;
};

class Database {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & carboncopy;
		ar & names;
		ar & masses;
		ar & volumes;
		ar & sfactors;
		ar & atomIDs;
    }
	/////////////////// 
	Database() {}
	Database(const Database&);
	Database& operator=(const Database&);

	std::vector<std::string> carboncopy;

	void read_conf(std::string filename);

	std::string guessformat(std::string filename);
	bool check();

public:
	DatabaseNamesParameters names;
	DatabaseMassesParameters masses;
	DatabaseVolumesParameters volumes;	
	DatabaseSFactorsParameters sfactors;	
	DatabaseAtomIDsParameters atomIDs;	

	void init(std::string filename, std::string format="");
	
	// interface for initiatilzation and interfacing
	static Database* Inst() { static Database instance; return &instance;}
	
	~Database() {}; // it is said some compilers have problems w/ private destructors.
	
	void write(std::string filename, std::string format="");
};

#endif 