/*
 *  coordinate_sets.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */
// direct header
#include "sample/coordinate_sets.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "sample/frame.hpp"
#include "coor3d.hpp"
#include "measures/center_of_mass.hpp"
#include "control.hpp"


using namespace std;

// has to be default constructible
CoordinateSets::CoordinateSets() {

}

CoordinateSets::~CoordinateSets() {
	typedef std::map<size_t,CoordinateSet*>::iterator sci_iterator;
	for(sci_iterator sci = setcache.begin(); sci != setcache.end(); ++sci)
	{
		delete sci->second;
	}
	
	for(size_t i = 0; i < m_motion_walkers.size(); ++i)
	{
		delete m_motion_walkers[i].second;
	}
}

void CoordinateSets::init() {


    // read in frame information
    for(size_t i = 0; i < Params::Inst()->sample.frames.size(); ++i)
    {
    	SampleFramesetParameters& f = Params::Inst()->sample.frames[i];
    	Info::Inst()->write(string("Reading frames from: ")+f.filename);
    	size_t nof = frames.add_frameset(f.filename,f.type,f.first,f.last,f.last_set,f.stride);			
    	Info::Inst()->write(string("Found ")+to_s(nof)+string(" frames"));			
    }
    
	// construct a helper for the motions:
	if (Params::Inst()->sample.motions.size()>0) {
		for(size_t i = 0; i < Params::Inst()->sample.motions.size(); ++i)
		{
			SampleMotionParameters& motion = Params::Inst()->sample.motions[i];
			MotionWalker* p_mw = NULL;
			if (motion.type=="linear") {
				p_mw = new LinearMotionWalker(motion.displace,motion.direction);
			}
			else if (motion.type=="fixed") {
				p_mw = new FixedMotionWalker(motion.displace,motion.direction);
			}
			else if (motion.type=="oscillation") {
				p_mw = new OscillationMotionWalker(motion.displace,motion.frequency, motion.direction);
			}
			else if (motion.type=="randomwalk") {
				p_mw = new RandomMotionWalker(motion.displace,motion.seed, motion.direction);				
			}
			else if (motion.type=="brownian") {
				p_mw = new BrownianMotionWalker(motion.displace,motion.seed, motion.direction);				
			}
			else if (motion.type=="none") {
                p_mw = NULL;
			}
			else {
                Err::Inst()->write("Motion type not understood");
                throw;
			}
			
			
			if (p_mw!=NULL) {
				m_motion_walkers.push_back(make_pair<string,MotionWalker*>(motion.selection,p_mw));				
			}
		}
	}
	
    
	// construct a helper for the alignments:
	if (Params::Inst()->sample.alignments.size()>0) {
		for(size_t i = 0; i < Params::Inst()->sample.alignments.size(); ++i)
		{
            std::string type = Params::Inst()->sample.alignments[i].type;
            std::string selection = Params::Inst()->sample.alignments[i].selection;
            
            if (Params::Inst()->sample.alignments[i].order=="pre") {
                m_prealignments.push_back(make_pair(selection,type));
            } else if (Params::Inst()->sample.alignments[i].order=="post") {
                m_postalignments.push_back(make_pair(selection,type));                
            } else {
                Err::Inst()->write("Ordering of alignment not understood. Must be pre or post.");
                throw;
            }
            
		}
	}	
    
}
void CoordinateSets::load_into_cache(size_t framenumber) {
	CartesianCoordinateSet* pcset = NULL;
	if ( Params::Inst()->runtime.limits.cache.coordinate_sets < setcache.size() ) {
		clear_cache();
	}
	Frame frame = frames.load(framenumber);
	pcset = new CartesianCoordinateSet(frame,*p_selection);

    // align here
    for(size_t i = 0; i < m_prealignments.size(); ++i)
    {
		string sel = m_prealignments[i].first;
		string type = m_prealignments[i].second;
        if (type=="center") {
			p_atoms->assert_selection(sel);        
            CartesianCoor3D origin = CenterOfMass(*p_atoms,p_atoms->selections[sel],*p_selection,*pcset);
            pcset->translate(-1.0*origin,*p_selection, p_atoms->selections[sel]);
        }
    }
    
	// apply motions, operate in cartesian space
	for(size_t i = 0; i < m_motion_walkers.size(); ++i)
	{
		string sel = m_motion_walkers[i].first;
		MotionWalker* p_mw = m_motion_walkers[i].second;

		if (sel == "") {
			pcset->translate(p_mw->translation(framenumber),*p_selection, p_atoms->selections[sel]);		
		} else {
			p_atoms->assert_selection(sel);
			pcset->translate(p_mw->translation(framenumber),*p_selection, p_atoms->selections[sel]);						
		}	
	}		

    // align here
    for(size_t i = 0; i < m_postalignments.size(); ++i)
    {
		string sel = m_postalignments[i].first;
		string type = m_postalignments[i].second;
        if (type=="center") {
			p_atoms->assert_selection(sel);        
            CartesianCoor3D origin = CenterOfMass(*p_atoms,p_atoms->selections[sel],*p_selection,*pcset);
            pcset->translate(-1.0*origin,*p_selection, p_atoms->selections[sel]);
        }
    }
    
	if (m_representation==CARTESIAN) {
 	    setcache[framenumber] = pcset;
	} else if (m_representation==SPHERICAL) {
	    setcache[framenumber] = new SphericalCoordinateSet(*pcset); 	    
        delete pcset;
	} else if (m_representation==CYLINDRICAL) {
	    setcache[framenumber] = new CylindricalCoordinateSet(*pcset); 	    
        delete pcset;
	}	
}

void CoordinateSets::set_representation(CoordinateRepresentation representation) {
    // a change of representation invalides the cache:
    if (m_representation!=representation) {
        clear_cache();        
    }
    
    m_representation = representation;
}

CoordinateRepresentation CoordinateSets::get_representation() {
    return m_representation;
}

void CoordinateSets::set_atoms(Atoms& atoms) {
    p_atoms = &atoms;
}

CoordinateSet& CoordinateSets::load(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	// this is the cache:
	if (setcache.find(framenumber)==setcache.end()) {	
		load_into_cache(framenumber);
	}

	pcset = setcache[framenumber];
	
	currentframe_i = framenumber;
	return *pcset;
}

void CoordinateSets::write_xyz(std::string filename) {
    ofstream ofile(filename.c_str());
	for(size_t i = 0; i < setcache.size(); ++i) {
        CoordinateSet* pcset = setcache[i];
        ofile << pcset->size() << endl;
        ofile << "generated by s_coordump" << endl;
        for(size_t j = 0; j < pcset->size(); ++j)
        {
            ofile << j << " " << pcset->c1[j] << " " << pcset->c2[j] << " " << pcset->c3[j] << endl;
        }
	}
    ofile.close();
}

void CoordinateSets::clear_cache() {
	typedef std::map<size_t,CoordinateSet*>::iterator sci_iterator;
	for(sci_iterator sci = setcache.begin(); sci != setcache.end(); ++sci)
	{
		delete sci->second;
	}
	setcache.clear();
}

void CoordinateSets::set_selection(Atomselection& selection) {
    if (p_selection!=&selection) {
    	clear_cache();        
    }
	p_selection = &selection;
}

Atomselection& CoordinateSets::get_selection() {
	return *p_selection;
}

CoordinateSet& CoordinateSets::current() {
	return *(setcache[currentframe_i]);
}



// end of file
