/*
 *  coordinate_sets.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SAMPLE__COORDINATESETS_HPP_
#define SAMPLE__COORDINATESETS_HPP_

// common header
#include "common.hpp"

// standard header
#include <string>
#include <vector>

// special library headers

// other headers
#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "sample/frames.hpp"
#include "sample/coordinate_set.hpp"
#include "sample/motion_walker.hpp"

//forward declaration...

////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////

// This class is used by Scatterdevices to manage coordinate sets
class CoordinateSets {
protected:
    friend class boost::serialization::access;	
	template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
        ar.register_type(static_cast<CoordinateSet*>(NULL));
        ar.register_type(static_cast<Atoms*>(NULL));
        ar.register_type(static_cast<Atomselection*>(NULL));
        
        ar.register_type(static_cast<LinearMotionWalker*>(NULL));
        ar.register_type(static_cast<RandomMotionWalker*>(NULL));
        ar.register_type(static_cast<FixedMotionWalker*>(NULL));
        ar.register_type(static_cast<OscillationMotionWalker*>(NULL));
        ar.register_type(static_cast<BrownianMotionWalker*>(NULL));
        
		ar & setcache;
		ar & currentframe_i;
		ar & frames;
        ar & m_motion_walkers;
        ar & m_prealignments;
        ar & m_postalignments;
        ar & m_representation;
        
		// DONT set back references. have to be set from the outside
		// p_atoms;
		// p_selection;
    }
    
	std::map<size_t,CoordinateSet*> setcache;
	std::vector< std::pair<std::string,MotionWalker*> > m_motion_walkers;
	std::vector< std::pair<std::string,std::string> > m_prealignments;
	std::vector< std::pair<std::string,std::string> > m_postalignments;

    Frames frames;

	Atoms* p_atoms;
	Atomselection* p_selection;
	
	size_t currentframe_i;

    CoordinateRepresentation m_representation;
	void load_into_cache(size_t framenumber);
	
public:
	CoordinateSets();
	~CoordinateSets() ;
	
	CoordinateSet& current();
	CoordinateSet& load(size_t frame);	
	void clear_cache();
	
	// use these to initialize the coordinate set:
	void set_selection(Atomselection& selection);
	void set_atoms(Atoms& atoms);
	
    void set_representation(CoordinateRepresentation representation);
    CoordinateRepresentation get_representation();
	
    void write_xyz(std::string filename); // dumps coordinates to file in xyz format
	
	Atomselection& get_selection();
    
    void init();
    size_t size() { return frames.size(); }
};

#endif
// end of file