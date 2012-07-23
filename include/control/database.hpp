/** \file
This file contains the interface for the application wide database, which, among other things, defines physical constants and atom name mapping.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef CONTROL__DATABASE_HPP_
#define CONTROL__DATABASE_HPP_

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

// This class follow a strong hierarchy. The root class is located at the end of this file.

/**
Implements a mapping between integer based atom IDs and exclusion factors. Exclusion factors are a means to scale atomic volumes.
*/
class DatabaseExlusionParameters {
friend class Database;      
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
	DatabaseExlusionParameters() : quicklookup_counter(0) {}

//	void reg(std::string label, std::vector<double> constants,size_t function_type);
//	double get(std::string label);
	void reg(size_t ID, std::vector<double> constants,size_t function_type);
	double get(size_t ID,double effvolume,double q);

};

/**
Implements a mapping between integer based atom IDs and their associated volume.
*/
class DatabaseVolumesParameters {
friend class Database;         
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

/**
Implements a mapping between integer based atom IDs and their potentially q dependent scattering factors.
*/
class DatabaseSFactorsParameters {
friend class Database;     
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

/**
Implements a mapping between integer based atom IDs and their masses.
*/
class DatabaseMassesParameters {
friend class Database;         
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

/**
Implements a mapping between integer based atom IDs and string labels. For efficency, internally all physical constants are associated with integer based atom IDs.
*/
class DatabaseAtomIDsParameters {
friend class Database;    
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

/**
Support for PDB style atom names
*/
class DatabaseNamesPDBParameters {
friend class Database;        
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

/**
Implements the name mapping facility of the database
*/
class DatabaseNamesParameters {
friend class Database;     
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

/**
This singleton class is the root of a strong class hierarchy which implements a software wide database for physical constants and atom labeling. It also performs required computation, e.g. when the scattering factor for an atom is q vector length dependent. 
*/
class Database {
private:
	/////////////////// MPI related
	// make this class serializable to 
	// allow sample to be transmitted via MPI
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
		ar & rawconfig;
		ar & config;
		ar & names;
		ar & masses;
		ar & volumes;
		ar & exclusionfactors;
		ar & sfactors;
		ar & atomIDs;
    }
	/////////////////// 
	Database() {}
	Database(const Database&);
	Database& operator=(const Database&);

	std::vector<char> rawconfig;
	std::vector<char> config;

	void read_xml(std::string filename);

	std::string guessformat(std::string filename);
	bool check();

public:
	DatabaseNamesParameters names;
	DatabaseMassesParameters masses;
	DatabaseVolumesParameters volumes;	
	DatabaseExlusionParameters exclusionfactors;	
	DatabaseSFactorsParameters sfactors;	
	DatabaseAtomIDsParameters atomIDs;	

	void init();
	
	void get_rawconfig(std::vector<char>& rc) { rc = rawconfig; }
	void get_config(std::vector<char>& c) { c=config; }
		
	// interface for initiatilzation and interfacing
	static Database* Inst() { static Database instance; return &instance;}
	
	~Database() {}; // it is said some compilers have problems w/ private destructors.
	
	void write(std::string filename, std::string format="");
};

#endif

// end of file
