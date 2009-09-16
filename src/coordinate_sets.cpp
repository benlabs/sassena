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
#include "parameters.hpp"
#include "database.hpp"

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
}

CoordinateSet& CoordinateSets::load(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	if (setcache.find(framenumber)==setcache.end()) {
		if ((Params::Inst()->limits.coordinatesets_cache_max>0) && setcache.size()>=Params::Inst()->limits.coordinatesets_cache_max) {
			clear_cache();
		}
		
		p_sample->frames.load(framenumber,p_sample->atoms);
		pcset = new CoordinateSet(p_sample->frames.current(),*p_selection);
		setcache[framenumber] = pcset;
		
		if (Params::Inst()->sample.motions.size()>0) {
			for(size_t i = 0; i < Params::Inst()->sample.motions.size(); ++i)
			{
				SampleMotionParameters& motion = Params::Inst()->sample.motions[i];
				if (motion.type=="linear") {
					pcset->translate(framenumber*motion.displace*motion.direction);
				}
				else if (motion.type=="fixed") {
					pcset->translate(motion.displace*motion.direction);
				}
				else if (motion.type=="oscillation") {
					pcset->translate(motion.displace*motion.direction*sin(framenumber*motion.frequency));
				}
				
			}
		}
	}
	else {
		pcset = setcache[framenumber];
	}
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




// has to be default constructible
CoordinateSetsM::CoordinateSetsM() {
	p_origin_selection = NULL;		
}

CoordinateSetsM::~CoordinateSetsM() {
	typedef std::map<size_t,CoordinateSetM*>::iterator sci_iterator;
	for(sci_iterator sci = setcache.begin(); sci != setcache.end(); ++sci)
	{
		delete sci->second;
	}
}

CoordinateSetM& CoordinateSetsM::load(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	CoordinateSetM* pcsetm = NULL;
	
	if (setcache.find(framenumber)==setcache.end()) {
		if ((Params::Inst()->limits.coordinatesets_cache_max>0) && setcache.size()>=Params::Inst()->limits.coordinatesets_cache_max) {
			clear_cache();
		}
		
		p_sample->frames.load(framenumber,p_sample->atoms);
		pcset = new CoordinateSet(p_sample->frames.current(),*p_selection);
		
		if (Params::Inst()->sample.motions.size()>0) {
			for(size_t i = 0; i < Params::Inst()->sample.motions.size(); ++i)
			{
				SampleMotionParameters& motion = Params::Inst()->sample.motions[i];
				if (motion.type=="linear") {
					pcset->translate(framenumber*motion.displace*motion.direction);
				}
				else if (motion.type=="fixed") {
					pcset->translate(motion.displace*motion.direction);
				}
				else if (motion.type=="oscillation") {
					pcset->translate(motion.displace*motion.direction*sin(framenumber*motion.frequency));
				}
				
			}
		}

		CartesianCoor3D origin(0,0,0);
		if (p_origin_selection!=NULL) {
			origin = p_sample->frames.current().cofm(p_sample->atoms,*p_origin_selection);
		}
		
		pcsetm = new CoordinateSetM(pcset,origin);
		delete pcset;

		setcache[framenumber] = pcsetm;
	}
	else {
		pcsetm = setcache[framenumber];
	}
	currentframe_i = framenumber;
	return *pcsetm;
}

void CoordinateSetsM::clear_cache() {
	typedef std::map<size_t,CoordinateSetM*>::iterator sci_iterator;
	for(sci_iterator sci = setcache.begin(); sci != setcache.end(); ++sci)
	{
		delete sci->second;
	}
	setcache.clear();
}

void CoordinateSetsM::set_selection(Atomselection& selection) {
	clear_cache();
	p_selection = &selection;
}

Atomselection& CoordinateSetsM::get_selection() {
	return *p_selection;
}

void CoordinateSetsM::set_sample(Sample& sample) {
	clear_cache();
	p_sample = &sample;
}

void CoordinateSetsM::set_origin(Atomselection& origin) {
	p_origin_selection = &origin;
}

CoordinateSetM& CoordinateSetsM::current() {
	return *(setcache[currentframe_i]);
}


// end of file