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
#include <sstream>

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
		ar & filepath;
		ar & format;
    }
	/////////////////// 

public:	
	std::string file;
	std::string filepath; // runtime
	std::string format;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;
		ss << std::string(pad,' ') << "<format>" << format << "</format>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type_ << "</type>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << SampleSelectionParameters::write_xml(pad);
		for(size_t i = 0; i < ids_.size(); ++i)
		{
			ss << std::string(pad,' ') << "<index>" << ids_[i] << "</index>" << std::endl;	
		}
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << SampleSelectionParameters::write_xml(pad);
		ss << std::string(pad,' ') << "<from>" << from_ << "</from>" << std::endl;	
		ss << std::string(pad,' ') << "<to>" << to_ << "</to>" << std::endl;	
		return ss.str();
	}

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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << SampleSelectionParameters::write_xml(pad);
		ss << std::string(pad,' ') << "<expression>" << expression_ << "</expression>" << std::endl;	
		return ss.str();
	}
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
	    ar & filepath_;
        ar & format_;
        ar & selector_;
        ar & expression_;
    }
public:

    std::string file_;
    std::string filepath_; // runtime
    std::string format_;
    std::string selector_;    
    std::string expression_;
    
    SampleFileSelectionParameters() {}
    SampleFileSelectionParameters(std::string file, std::string format, std::string selector, std::string expression) : SampleSelectionParameters("file"), file_(file), format_(format), selector_(selector), expression_(expression) {}
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << SampleSelectionParameters::write_xml(pad);
		ss << std::string(pad,' ') << "<file>" << file_ << "</file>" << std::endl;	
		ss << std::string(pad,' ') << "<format>" << format_ << "</format>" << std::endl;	
		ss << std::string(pad,' ') << "<selector>" << selector_ << "</selector>" << std::endl;	
		ss << std::string(pad,' ') << "<expression>" << expression_ << "</expression>" << std::endl;	
		return ss.str();
	}
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
		ar & file;
		ar & filepath;
		ar & format;
        ar & clones;
        ar & index;
        ar & indexpath;
		ar & index_default;
    }
	/////////////////// 

public:	
	size_t first;
	size_t last;
    size_t clones;
	bool last_set;
	size_t stride;
	std::string file;
	std::string filepath; // runtime
	std::string format;
    std::string index;
    std::string indexpath;
	bool index_default;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<format>" << format << "</format>" << std::endl;	
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;	
		// don't write location specific filepath or indexpath! not well defined behavior
		if (!index_default) {
			ss << std::string(pad,' ') << "<index>" << index << "</index>" << std::endl;	
		}
		ss << std::string(pad,' ') << "<first>" << first << "</first>" << std::endl;
		if (last_set) {
			ss << std::string(pad,' ') << "<last>" << last << "</last>" << std::endl;				
		}	
		ss << std::string(pad,' ') << "<clones>" << clones << "</clones>" << std::endl;	
		ss << std::string(pad,' ') << "<stride>" << stride << "</stride>" << std::endl;	
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		for(size_t i = 0; i < this->size(); ++i)
		{
			ss << std::string(pad,' ') << "<frameset>" << std::endl;
			ss << this->at(i).write_xml(pad+1);
			ss << std::string(pad,' ') << "</frameset>" << std::endl;	
		}
		return ss.str();
	}
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
        ar & filepath;
        ar & format;
        ar & selection;
    }
	/////////////////// 

public:	
	std::string type;
	std::string selection;
	std::string file;
	std::string filepath; // runtime
	std::string format;
	size_t frame;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<selection>" << selection << "</selection>" << std::endl;	
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;	
		ss << std::string(pad,' ') << "<format>" << format << "</format>" << std::endl;	
		ss << std::string(pad,' ') << "<frame>" << frame << "</frame>" << std::endl;	
		return ss.str();
	}
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
	unsigned long seed;
	long sampling;	
	CartesianCoor3D direction;
    SampleMotionReferenceParameters reference;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<displace>" << displace << "</displace>" << std::endl;	
		ss << std::string(pad,' ') << "<frequency>" << frequency << "</frequency>" << std::endl;	
		ss << std::string(pad,' ') << "<radius>" << radius << "</radius>" << std::endl;	
		ss << std::string(pad,' ') << "<selection>" << selection << "</selection>" << std::endl;	
		ss << std::string(pad,' ') << "<seed>" << seed << "</seed>" << std::endl;	
		ss << std::string(pad,' ') << "<sampling>" << sampling << "</sampling>" << std::endl;	
		ss << std::string(pad,' ') << "<referemce>" << std::endl;
		ss << reference.write_xml(pad+1);
		ss << std::string(pad,' ') << "</referemce>" << std::endl;
		return ss.str();
	}
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
        ar & filepath;
        ar & format;
        ar & selection;
    }
	/////////////////// 

public:	
	std::string type;
	std::string selection;
	std::string file;
	std::string filepath; // runtime	
	std::string format;
	size_t frame;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<selection>" << selection << "</selection>" << std::endl;	
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;	
		ss << std::string(pad,' ') << "<format>" << format << "</format>" << std::endl;	
		ss << std::string(pad,' ') << "<frame>" << frame << "</frame>" << std::endl;	
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<selection>" << selection << "</selection>" << std::endl;	
		ss << std::string(pad,' ') << "<order>" << order << "</order>" << std::endl;	
		ss << std::string(pad,' ') << "<referemce>" << std::endl;
		ss << reference.write_xml(pad+1);
		ss << std::string(pad,' ') << "</referemce>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<structure>" << std::endl;
		ss << structure.write_xml(pad+1);
		ss << std::string(pad,' ') << "</structure>" << std::endl;
		ss << std::string(pad,' ') << "<framesets>" << std::endl;
		ss << framesets.write_xml(pad+1);
		ss << std::string(pad,' ') << "</framesets>" << std::endl;		
		ss << std::string(pad,' ') << "<selections>" << std::endl;
		for(std::map<std::string,SampleSelectionParameters*>::iterator i = selections.begin(); i != selections.end(); ++i)
		{
			ss << std::string(pad+1,' ') << "<selection>" << std::endl;
			ss << std::string(pad+2,' ') << "<name>" << i->first << "</name>" << std::endl;
			ss << i->second->write_xml(pad+2);
			ss << std::string(pad+1,' ') << "</selection>" << std::endl;			
		}
		ss << std::string(pad,' ') << "</selections>" << std::endl;
		ss << std::string(pad,' ') << "<motions>" << std::endl;
		for(size_t i = 0; i < motions.size(); ++i)
		{
			ss << std::string(pad+1,' ') << "<motion>" << std::endl;
			ss << motions[i].write_xml(pad+2);
			ss << std::string(pad+1,' ') << "</motion>" << std::endl;			
		}
		ss << std::string(pad,' ') << "</motions>" << std::endl;
		ss << std::string(pad,' ') << "<alignments>" << std::endl;
		for(size_t i = 0; i < alignments.size(); ++i)
		{
			ss << std::string(pad+1,' ') << "<alignment>" << std::endl;
			ss << alignments[i].write_xml(pad+2);
			ss << std::string(pad+1,' ') << "</alignment>" << std::endl;			
		}
		ss << std::string(pad,' ') << "</alignments>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<selection>" << selection << "</selection>" << std::endl;	
		ss << std::string(pad,' ') << "<value>" << value << "</value>" << std::endl;	
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<factor>" << factor << "</factor>" << std::endl;	
		ss << std::string(pad,' ') << "<kappas>" << std::endl;	
		for(size_t i = 0; i < kappas.size(); ++i)
		{
			ss << std::string(pad+1,' ') << "<kappa>" << std::endl;	
			ss << kappas[i].write_xml(pad+2);
			ss << std::string(pad+1,' ') << "</kappa>" << std::endl;	
		}
		ss << std::string(pad,' ') << "</kappas>" << std::endl;	
		return ss.str();
	}
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
		ar & filepath;		
		ar & seed;
    }
	/////////////////// 

public:	
	std::string type;
	std::string algorithm;
	std::string file;
	std::string filepath; // runtime
	size_t resolution;
	unsigned long seed;
	
	void create();
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<algorithm>" << algorithm << "</algorithm>" << std::endl;	
		ss << std::string(pad,' ') << "<resolution>" << resolution << "</resolution>" << std::endl;	
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;	
		ss << std::string(pad,' ') << "<seed>" << seed << "</seed>" << std::endl;	
		return ss.str();
	}
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
        ar & filepath;
    }
	/////////////////// 

public:	
	std::string type;
    long resolution;
    std::string file;
    std::string filepath; // runtime

	void create();
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;	
		ss << std::string(pad,' ') << "<resolution>" << resolution << "</resolution>" << std::endl;	
		ss << std::string(pad,' ') << "<file>" << file << "</file>" << std::endl;	
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;
		
		ss << moments.write_xml(pad+1);

		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;
		return ss.str();
	}
	
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;
		ss << std::string(pad,' ') << "<axis>" << std::endl;
		ss << std::string(pad+1,' ') << "<x>" << axis.x << "</x>" << std::endl;
		ss << std::string(pad+1,' ') << "<y>" << axis.y << "</y>" << std::endl;
		ss << std::string(pad+1,' ') << "<z>" << axis.z << "</z>" << std::endl;
		ss << std::string(pad,' ') << "</axis>" << std::endl;
		ss << std::string(pad,' ') << "<vectors>" << std::endl;
		ss << vectors.write_xml(pad+1);
		ss << std::string(pad,' ') << "</vectors>" << std::endl;
		ss << std::string(pad,' ') << "<multipole>" << std::endl;
		ss << multipole.write_xml(pad+1);
		ss << std::string(pad,' ') << "</multipole>" << std::endl;
		ss << std::string(pad,' ') << "<exact>" << std::endl;
		ss << exact.write_xml(pad+1);
		ss << std::string(pad,' ') << "</exact>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<orientation>" << std::endl;
		ss << orientation.write_xml(pad+1);
		ss << std::string(pad,' ') << "</orientation>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>" << std::endl;
		ss << std::string(pad,' ') << "<method>" << method << "</method>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<from>" << from << "</from>" << std::endl;
		ss << std::string(pad,' ') << "<to>" << to << "</to>" << std::endl;
		ss << std::string(pad,' ') << "<points>" << points << "</points>" << std::endl;
		ss << std::string(pad,' ') << "<exponent>" << exponent << "</exponent>" << std::endl;
		ss << std::string(pad,' ') << "<base>" << std::endl;
		ss << std::string(pad,' ') << "<x>" << basevector.x << "</x>" << std::endl;
		ss << std::string(pad,' ') << "<y>" << basevector.y << "</y>" << std::endl;
		ss << std::string(pad,' ') << "<z>" << basevector.z << "</z>" << std::endl;
		ss << std::string(pad,' ') << "</base>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<scans>" << std::endl;
		for(size_t i = 0; i < scans.size(); ++i)
		{
			ss << std::string(pad+1,' ') << "<scan>" << std::endl;			
			ss << scans[i].write_xml(pad+2);
			ss << std::string(pad+1,' ') << "</scan>" << std::endl;			
		}
		ss << std::string(pad,' ') << "</scans>" << std::endl;
		return ss.str();
	}
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
		ar & filepath;		
        ar & fqt;
        ar & fq;
        ar & fq0;
        ar & fq2;
    }
	/////////////////// 

public:
	std::string file;
	std::string filepath; // runtime	
    bool fqt;
    bool fq;
    bool fq0;
    bool fq2;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<file>" << file << "</file>"<< std::endl;
		ss << std::string(pad,' ') << "<fqt>" << fqt << "</fqt>"<< std::endl;
		ss << std::string(pad,' ') << "<fq>" << fq << "</fq>"<< std::endl;
		ss << std::string(pad,' ') << "<fq0>" << fq0 << "</fq0>"<< std::endl;
		ss << std::string(pad,' ') << "<fq2>" << fq2 << "</fq2>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>"<< std::endl;
		ss << std::string(pad,' ') << "<dsp>" << std::endl;
		ss << dsp.write_xml(pad+1);
		ss << std::string(pad,' ') << "</dsp>" << std::endl;
		ss << std::string(pad,' ') << "<vectors>" << std::endl;
		ss << qvectors.write_xml(pad+1);
		ss << std::string(pad,' ') << "</vectors>" << std::endl;
		ss << std::string(pad,' ') << "<average>" << std::endl;
		ss << average.write_xml(pad+1);
		ss << std::string(pad,' ') << "</average>" << std::endl;
		ss << std::string(pad,' ') << "<background>" << std::endl;
		ss << background.write_xml(pad+1);
		ss << std::string(pad,' ') << "</background>" << std::endl;
		ss << std::string(pad,' ') << "<signal>" << std::endl;
		ss << signal.write_xml(pad+1);
		ss << std::string(pad,' ') << "</signal>" << std::endl;
		return ss.str();
	}
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
		ar & filepath;		
		ar & format;
		ar & target;
        ar & mode;
    }
	/////////////////// 

public:	
    bool dump;
	std::string file;
	std::string filepath; // runtime	
	std::string format;
	std::string target;
    std::string mode;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<dump>" << dump << "</dump>"<< std::endl;
		ss << std::string(pad,' ') << "<file>" << file << "</file>"<< std::endl;
		ss << std::string(pad,' ') << "<format>" << format << "</format>"<< std::endl;
		ss << std::string(pad,' ') << "<target>" << target << "</target>"<< std::endl;
		ss << std::string(pad,' ') << "<mode>" << mode << "</mode>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<server>" << server << "</server>"<< std::endl;
		ss << std::string(pad,' ') << "<client>" << client << "</client>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<serverflush>" << serverflush << "</serverflush>"<< std::endl;
		ss << std::string(pad,' ') << "<clientflush>" << clientflush << "</clientflush>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<memory>" << std::endl;
		ss << memory.write_xml(pad+1);		
		ss << std::string(pad,' ') << "</memory>" << std::endl;
		ss << std::string(pad,' ') << "<times>" << std::endl;
		ss << times.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</times>" << std::endl;						
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<delay>" << delay << "</delay>"<< std::endl;
		ss << std::string(pad,' ') << "<sampling>" << sampling << "</sampling>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<signal>" << std::endl;
		ss << signal.write_xml(pad+1);		
		ss << std::string(pad,' ') << "</signal>" << std::endl;
		ss << std::string(pad,' ') << "<monitor>" << std::endl;
		ss << monitor.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</monitor>" << std::endl;						
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<result_buffer>" << result_buffer << "</result_buffer>"<< std::endl;
		ss << std::string(pad,' ') << "<alignpad_buffer>" << alignpad_buffer << "</alignpad_buffer>"<< std::endl;
		ss << std::string(pad,' ') << "<exchange_buffer>" << exchange_buffer << "</exchange_buffer>"<< std::endl;
		ss << std::string(pad,' ') << "<signal_buffer>" << signal_buffer << "</signal_buffer>"<< std::endl;
		ss << std::string(pad,' ') << "<scale>" << scale << "</scale>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<threads>" << threads << "</threads>"<< std::endl;
		ss << std::string(pad,' ') << "<processes>" << processes << "</processes>"<< std::endl;
		ss << std::string(pad,' ') << "<cores>" << cores << "</cores>"<< std::endl;
		ss << std::string(pad,' ') << "<memory>" << std::endl;
		ss << memory.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</memory>" << std::endl;
		return ss.str();
	}
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
    std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<automatic>" << automatic << "</automatic>"<< std::endl;
		ss << std::string(pad,' ') << "<size>" << size << "</size>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<utilization>" << utilization << "</utilization>"<< std::endl;
		ss << std::string(pad,' ') << "<partitions>" << std::endl;
		ss << partitions.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</partitions>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<chunksize>" << chunksize << "</chunksize>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<data>" << data << "</data>"<< std::endl;
		ss << std::string(pad,' ') << "<buffer>" << buffer << "</buffer>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<memory>" << std::endl;
		ss << memory.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</memory>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<stage>" << std::endl;
		ss << stage.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</stage>" << std::endl;
		ss << std::string(pad,' ') << "<signal>" << std::endl;
		ss << signal.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</signal>" << std::endl;
		ss << std::string(pad,' ') << "<services>" << std::endl;
		ss << services.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</services>" << std::endl;
		ss << std::string(pad,' ') << "<computation>" << std::endl;
		ss << computation.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</computation>" << std::endl;
		ss << std::string(pad,' ') << "<decomposition>" << std::endl;
		ss << decomposition.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</decomposition>" << std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<update>" << update << "</update>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<orientations>" << orientations << "</orientations>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<write>" << write << "</write>"<< std::endl;
		ss << std::string(pad,' ') << "<buffer>" << buffer << "</buffer>"<< std::endl;
		return ss.str();
	}
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
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<timer>" << timer << "</timer>"<< std::endl;
		ss << std::string(pad,' ') << "<barriers>" << barriers << "</barriers>"<< std::endl;
		ss << std::string(pad,' ') << "<iowrite>" << std::endl;
		ss << iowrite.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</iowrite>" << std::endl;
		ss << std::string(pad,' ') << "<print>" << std::endl;
		ss << print.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</print>" << std::endl;
		ss << std::string(pad,' ') << "<monitor>" << std::endl;
		ss << monitor.write_xml(pad+1);	
		ss << std::string(pad,' ') << "</monitor>" << std::endl;
		return ss.str();
	}
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
		ar & filepath;		
        ar & format;
    }
	/////////////////// 

public:
	std::string type;
	std::string file;
	std::string filepath;	
	std::string format;
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<type>" << type << "</type>"<< std::endl;
		ss << std::string(pad,' ') << "<file>" << file << "</file>"<< std::endl;
		ss << std::string(pad,' ') << "<format>" << format << "</format>"<< std::endl;
		return ss.str();
	}
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
		ar & rawconfig;
		ar & config;
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

	std::vector<char> rawconfig; // contains raw copy of the input configuration
	std::vector<char> config; // contains complete text copy of the input configuration
	
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
	
	void get_rawconfig(std::vector<char>& rc) { rc = rawconfig; }
	void get_config(std::vector<char>& c) { c=config; }
	
	void init(int argc,char** argv);
	~Params() {}; // it is said some compilers have problems w/ private destructors.
	
	void write_xml_to_file(std::string filename);
	std::string write_xml(int pad=0) {
		std::stringstream ss;
		ss << std::string(pad,' ') << "<root>" << std::endl;
		
		ss << std::string(pad+1,' ') << "<sample>" << std::endl;
		ss << sample.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</sample>" << std::endl;
		ss << std::string(pad+1,' ') << "<scattering>" << std::endl;
		ss << scattering.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</scattering>" << std::endl;
		ss << std::string(pad+1,' ') << "<stager>" << std::endl;
		ss << stager.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</stager>" << std::endl;
		ss << std::string(pad+1,' ') << "<database>" << std::endl;
		ss << database.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</database>" << std::endl;
		ss << std::string(pad+1,' ') << "<limits>" << std::endl;
		ss << limits.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</limits>" << std::endl;
		ss << std::string(pad+1,' ') << "<debug>" << std::endl;
		ss << debug.write_xml(pad+2);	
		ss << std::string(pad+1,' ') << "</debug>" << std::endl;

		ss << std::string(pad,' ') << "</root>" << std::endl;
		return ss.str();
	}	
};

#endif 

// end of file
