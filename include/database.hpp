/*
 *  database.hpp
 *reg
 *  Created on: June 25, 2009
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef DATABASE_HPP_
#define DATABASE_HPP_

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
	
	size_t size() { return nextID; }
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
	void read_xml(std::string filename);

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

// end of file
