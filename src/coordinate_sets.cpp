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
	
	// apply alignment
	if (Params::Inst()->sample.alignment.origin.type=="manual") {
        pcset->translate(-1.0*Params::Inst()->sample.alignment.origin.basevector);
	} else if (Params::Inst()->sample.alignment.origin.type=="auto") {
        p_sample->atoms.assert_selection(Params::Inst()->sample.alignment.origin.selection);
        Atomselection& thissel = p_sample->atoms.selections[Params::Inst()->sample.alignment.origin.selection];
        CartesianCoor3D com = CenterOfMass(p_sample->atoms,p_sample->frames.current(),thissel);
		pcset->translate(-1.0*com);
	} else {
        Err::Inst()->write("sample.alignment.origin.type not understood");
        throw;
	}
	
	if (Params::Inst()->sample.alignment.axis.type=="manual") {
        pcset->rotate(Params::Inst()->sample.alignment.axis.basevector,CartesianCoor3D(0,1,0));
	} else if (Params::Inst()->sample.alignment.axis.type=="auto") {
        Err::Inst()->write("sample.alignment.axis.type==auto not implemented yet");
        throw;
	} else {
        Err::Inst()->write("sample.alignment.axis.type not understood");
        throw;
	}
	
	// apply motions
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

void CoordinateSets::write_xyz(std::string filename) {
    ofstream ofile(filename.c_str());
	for(size_t i = 0; i < setcache.size(); ++i) {
        CoordinateSet* pcset = setcache[i];
        ofile << pcset->size() << endl;
        ofile << "generate by s_coordump" << endl;
        for(size_t j = 0; j < pcset->size(); ++j)
        {
            ofile << j << " " << pcset->x[j] << " " << pcset->y[j] << " " << pcset->z[j] << endl;
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

CoordinateSet& CoordinateSetsMS::load(size_t framenumber) {
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

void CoordinateSetsMS::set_origin(Atomselection& origin) {
	if (! origin.is_subset_of(*p_selection)) {
		Err::Inst()->write("Origin selection must be a subset of CoordinateSet selection");
		throw;
	}
    p_origin_selection = &origin;
}


//// C:

CoordinateSet& CoordinateSetsMC::load(size_t framenumber) {
	CoordinateSet* pcset = NULL;
	// this is the cache:
	if (setcache.find(framenumber)==setcache.end()) {	
		CoordinateSets::load_into_cache(framenumber); // call parent function to handle class behaviour
		pcset = setcache[framenumber];

		CartesianCoor3D origin; origin = CenterOfMass(p_sample->atoms,*p_origin_selection,*p_selection,*pcset);
		// now do the transformation:
		CartesianCoor3D axis = Params::Inst()->scattering.average.orientation.axis;
		axis = axis / axis.length();
		
		// fix this:
        CartesianCoor3D qperpenticular;
        CartesianCoor3D o;
        double qperpenticular_l = 0;
		
		
		for(size_t i = 0; i < pcset->size(); ++i)
		{
			double x = pcset->x[i] - origin.x;
			double y = pcset->y[i] - origin.y;
			double z = pcset->z[i] - origin.z;
            CartesianCoor3D ct(x,y,z);
			
			CartesianCoor3D ctparallel = (axis*ct)*axis; 
			CartesianCoor3D ctperpenticular = ct - ctparallel; // this contains psi-phi
			double ctperpenticular_l = ctperpenticular.length(); // this is r
			double ctparallel_l = ctparallel.length(); // this is z
			
			double psiphi = 0;

			// if either qper_l or ctper_l is zero , then the bessel_terms vanish and delta_phipsi is irrelevant
			if ((qperpenticular_l!=0) && (ctperpenticular_l!=0)) {
				psiphi = acos( (ctperpenticular * qperpenticular) / (ctperpenticular_l*qperpenticular_l));
				CartesianCoor3D ctq = (ctperpenticular.cross_product(qperpenticular));
				ctq = ctq / ctq.length();				
				if ((ctq + o).length()>1.0) psiphi = 2*M_PI-psiphi; // if o || ctq -> 2; 0 otherwise			
			}
			
			CylinderCoor3D sc(CartesianCoor3D(x,y,z));
			// dirty hack: make x, y, z = r , phi , rho
			pcset->x[i] = sc.r;
			pcset->y[i] = sc.phi;
			pcset->z[i] = sc.z;
		}
	}

	pcset = setcache[framenumber];
	
	currentframe_i = framenumber;
	return *pcset;
}

void CoordinateSetsMC::set_origin(Atomselection& origin) {
	if (! origin.is_subset_of(*p_selection)) {
		Err::Inst()->write("Origin selection must be a subset of CoordinateSet selection");
		throw;
	}
    p_origin_selection = &origin;
}


// end of file