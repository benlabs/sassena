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
#include <boost/thread.hpp>

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
        
        ar.register_type(static_cast<LinearMotionWalker*>(NULL));
        ar.register_type(static_cast<RandomMotionWalker*>(NULL));
        ar.register_type(static_cast<FixedMotionWalker*>(NULL));
        ar.register_type(static_cast<OscillationMotionWalker*>(NULL));
        ar.register_type(static_cast<BrownianMotionWalker*>(NULL));
        ar.register_type(static_cast<LocalBrownianMotionWalker*>(NULL));
        
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

	std::vector< std::pair<std::string,MotionWalker*> > m_motion_walkers;
	std::vector< std::pair<std::string,std::string> > m_prealignments;
	std::vector< std::pair<std::string,std::string> > m_postalignments;

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
    void add_prealignment(std::string selection,std::string type);
    void add_postalignment(std::string selection,std::string type);
	
    void write_xyz(std::string filename); // dumps coordinates to file in xyz format
	
	IAtomselection& get_selection();
    
    void init();
    size_t size() { return frames.size(); }
};

#endif

// end of file
