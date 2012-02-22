/** \file
This file contains a class which defines a management class for coordinate sets.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef SAMPLE__COORDINATESETS_HPP_
#define SAMPLE__COORDINATESETS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers
#include <boost/thread.hpp>

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/frames.hpp"
#include "sample/coordinate_set.hpp"
#include "sample/motion_walker.hpp"

//forward declaration...

/** 
Stores alignment information which is used during the loading to perform transformation of the original trajectory data
*/
class CoordinateSetAlignment {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar & type;
        ar & order;    
        ar & selection; 
        ar & reference_selection; 
        ar & p_reference;
    }

public:
    std::string type;
    std::string order;    
    std::string selection; // atoms to move
    std::string reference_selection; // selection containing atoms used for fitting/alignment
    CoordinateSet* p_reference;
};

/** 
Stores motion walker alignment information which is used during the loading to combine the original trajectory data with artificial motions
*/
class MotionWalkerAlignment {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
	    ar.register_type(static_cast<LinearMotionWalker*>(NULL));
        ar.register_type(static_cast<RandomMotionWalker*>(NULL));
        ar.register_type(static_cast<FixedMotionWalker*>(NULL));
        ar.register_type(static_cast<OscillationMotionWalker*>(NULL));
        ar.register_type(static_cast<BrownianMotionWalker*>(NULL));
        ar.register_type(static_cast<LocalBrownianMotionWalker*>(NULL));
        ar.register_type(static_cast<RotationalBrownianMotionWalker*>(NULL));
    
        ar & p_mw;
        ar & selection; 
        ar & reference_selection; 
        ar & p_reference;
    }

public:
    MotionWalker* p_mw;
    std::string selection; // atoms to move
    std::string reference_selection; // selection containing atoms used for fitting/alignment
    CoordinateSet* p_reference;

};


/** 
Management class for coordinate sets.
*/
class CoordinateSets {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<CoordinateSet*>(NULL));
        ar.register_type(static_cast<Atoms*>(NULL));
                
        // ar & m_prealignmentvectors;
        // ar & m_postalignmentvectors;

        ar & m_motion_walkers;
        ar & m_prealignments;
        ar & m_postalignments;

		ar & frames;

		// DONT set back references. have to be set from the outside
		// p_atoms;
		// p_selection;

        ar & m_representation;
    }
    
    // std::map<size_t,std::vector<CartesianCoor3D> > m_prealignmentvectors;
    // std::map<size_t,std::vector<CartesianCoor3D> > m_postalignmentvectors;

	std::vector< MotionWalkerAlignment > m_motion_walkers;
	std::vector< CoordinateSetAlignment > m_prealignments;
	std::vector< CoordinateSetAlignment > m_postalignments;

    Frames frames;

	Atoms* p_atoms;
	IAtomselection* p_selection;
	
    CoordinateRepresentation m_representation;

public:
	CoordinateSets();
	~CoordinateSets() ;
	
	CoordinateSet* load(size_t frame);	
     
	// std::vector<CartesianCoor3D> get_prealignmentvectors(size_t frame);	
	// std::vector<CartesianCoor3D> get_postalignmentvectors(size_t frame);	
	
	// use these to initialize the coordinate set:
	void set_selection(IAtomselection* selection);
	void set_atoms(Atoms& atoms);
	
    void set_representation(CoordinateRepresentation representation);
    CoordinateRepresentation get_representation();
    void add_prealignment(CoordinateSetAlignment alignment);
    void add_postalignment(CoordinateSetAlignment alignment);
	
    void write_xyz(std::string filename); // dumps coordinates to file in xyz format
	
	IAtomselection& get_selection();
    
    void init();
    size_t size() { return frames.size(); }
};

#endif

// end of file
