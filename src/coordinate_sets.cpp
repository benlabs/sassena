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
					pcset->translate(motion.displace*motion.direction);
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

// end of file