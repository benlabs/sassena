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
#include "coordinate_sets.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// other headers
#include "frame.hpp"
#include "coor3d.hpp"
#include "log.hpp"
#include "measures/center_of_mass.hpp"
#include "parameters.hpp"
#include "database.hpp"

using namespace std;

// has to be default constructible
CoordinateSets::CoordinateSets() {
	
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
			
			if (p_mw!=NULL) {
				m_motion_walkers.push_back(make_pair<string,MotionWalker*>(motion.selection,p_mw));				
			}
		}
	}
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

void CoordinateSets::load_into_cache(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	if ((Params::Inst()->limits.coordinatesets_cache_max>0) && setcache.size()>=Params::Inst()->limits.coordinatesets_cache_max) {
		clear_cache();
	}
	p_sample->frames.load(framenumber,p_sample->atoms);
	pcset = new CoordinateSet(p_sample->frames.current(),*p_selection); 
	
	for(size_t i = 0; i < m_motion_walkers.size(); ++i)
	{
		string sel = m_motion_walkers[i].first;
		MotionWalker* p_mw = m_motion_walkers[i].second;

		if (sel == "") {
			pcset->translate(p_mw->translation(framenumber),*p_selection, p_sample->atoms.selections[sel]);		
		} else {
			p_sample->atoms.assert_selection(sel);
			pcset->translate(p_mw->translation(framenumber),*p_selection, p_sample->atoms.selections[sel]);						
		}	
	}		
	
	setcache[framenumber] = pcset;
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

void CoordinateSets::clear_cache() {
	typedef std::map<size_t,CoordinateSet*>::iterator sci_iterator;
	for(sci_iterator sci = setcache.begin(); sci != setcache.end(); ++sci)
	{
		delete sci->second;
	}
	setcache.clear();
}

void CoordinateSets::set_selection(Atomselection& selection) {
	clear_cache();
	p_selection = &selection;
}

Atomselection& CoordinateSets::get_selection() {
	return *p_selection;
}

void CoordinateSets::set_sample(Sample& sample) {
	clear_cache();
	p_sample = &sample;
}

CoordinateSet& CoordinateSets::current() {
	return *(setcache[currentframe_i]);
}


//// M:

CoordinateSet& CoordinateSetsM::load(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	// this is the cache:
	if (setcache.find(framenumber)==setcache.end()) {	
		CoordinateSets::load_into_cache(framenumber); // call parent function to handle class behaviour
		pcset = setcache[framenumber];

		CartesianCoor3D origin; origin = CenterOfMass(p_sample->atoms,*p_origin_selection,*p_selection,*pcset);
		// now do the transformation:
		for(size_t i = 0; i < pcset->size(); ++i)
		{
			double x = pcset->x[i] - origin.x;
			double y = pcset->y[i] - origin.y;
			double z = pcset->z[i] - origin.z;
			
			SphericalCoor3D sc(CartesianCoor3D(x,y,z));
			// dirty hack: make x, y, z = r , phi , rho
			pcset->x[i] = sc.r;
			pcset->y[i] = sc.phi;
			pcset->z[i] = sc.theta;
		}
	}

	pcset = setcache[framenumber];
	
	currentframe_i = framenumber;
	return *pcset;
}

void CoordinateSetsM::set_origin(Atomselection& origin) {
	if (! origin.is_subset_of(*p_selection)) {
		Err::Inst()->write("Origin selection must be a subset of CoordinateSet selection");
		throw;
	}
    p_origin_selection = &origin;
}

// end of file