/** \file
This file contains a class which defines a management class for coordinate sets.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "sample/coordinate_sets.hpp"

// standard header
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

// special library headers
#include <boost/lexical_cast.hpp>

// other headers
#include "sample/frame.hpp"
#include "math/coor3d.hpp"
#include "sample/center_of_mass.hpp"
#include "control.hpp"
#include "log.hpp"


using namespace std;

// has to be default constructible
CoordinateSets::CoordinateSets() {

}

CoordinateSets::~CoordinateSets() {
	for(size_t i = 0; i < m_motion_walkers.size(); ++i)
	{
		delete m_motion_walkers[i].p_reference;
		delete m_motion_walkers[i].p_mw;
	}
	for(size_t i = 0; i < m_prealignments.size(); ++i)
	{
        if (m_prealignments[i].p_reference!=NULL) {
            delete m_prealignments[i].p_reference;
        }
	}
	for(size_t i = 0; i < m_postalignments.size(); ++i)
	{
        if (m_postalignments[i].p_reference!=NULL) {
            delete m_postalignments[i].p_reference;
        }
	}
	
}

void CoordinateSets::init() {

    // read in frame information
    for(size_t i = 0; i < Params::Inst()->sample.framesets.size(); ++i)
    {
    	SampleFramesetParameters& f = Params::Inst()->sample.framesets[i];
    	Info::Inst()->write(string("Reading frames from: ")+f.filepath);
    	if (f.clones!=1) {
    	    Info::Inst()->write(string("Cloning frameset ")+boost::lexical_cast<string>(f.clones)+string(" times"));
    	}
        size_t nof = 0;
	    nof += frames.add_frameset(f.filepath,f.format,f.first,f.last,f.last_set,f.stride,f.indexpath,f.index_default,f.clones);			
    	Info::Inst()->write(string("Found ")+boost::lexical_cast<string>(nof)+string(" frames"));			
    }
	
	// construct a helper for the motion walker alignments:
	if (Params::Inst()->sample.motions.size()>0) {
		for(size_t i = 0; i < Params::Inst()->sample.motions.size(); ++i)
		{
			SampleMotionParameters& motion = Params::Inst()->sample.motions[i];
			MotionWalkerAlignment mwa;
			
			mwa.selection = motion.selection;
			mwa.reference_selection = motion.reference.selection;
			mwa.p_reference = NULL;
			mwa.p_mw = NULL;
			
            // read and store reference now
            if (motion.reference.type=="frame") {
                size_t framenumber = motion.reference.frame;
                if (framenumber>=this->size()) {
                    Warn::Inst()->write("reference frame number in alignment larger than size of frameset. Setting to last frame!");
                    framenumber = this->size()-1;
                }
                Frame frame = frames.load(framenumber);
                mwa.p_reference = new CartesianCoordinateSet(frame,p_atoms->selections[motion.reference.selection]);
            } else if (motion.reference.type=="file") {
                std::string file = motion.reference.file;
                std::string format = motion.reference.format;
                size_t framenumber = motion.reference.frame;
                if (format=="pdb") {
                    PDBFrameset fs(file,0);
                    fs.generate_index();
                    if (framenumber>=fs.number_of_frames) {
                        Warn::Inst()->write("reference frame number in alignment larger than size of frameset. Setting to last frame!");
                        framenumber = fs.number_of_frames-1;
                    }
                    Frame cf;		
                	fs.read_frame(framenumber,cf);
                	mwa.p_reference = new CartesianCoordinateSet(cf,p_atoms->selections[motion.reference.selection]);
                } else {
                    Err::Inst()->write(string("File format for alignment reference not understood, format=") + format);
                    throw;
                }
            } else if (motion.reference.type=="instant"){
				mwa.p_reference = NULL; // instructs code to use current coordinates
			} else {
                Err::Inst()->write(string("Reference type not understood:")+motion.reference.type );
                throw;
            }

			if (motion.type=="linear") {
				mwa.p_mw = new LinearMotionWalker(motion.displace,motion.sampling,motion.direction);
			}
			else if (motion.type=="fixed") {
				mwa.p_mw = new FixedMotionWalker(motion.displace,motion.direction);
			}
			else if (motion.type=="oscillation") {
				mwa.p_mw = new OscillationMotionWalker(motion.displace,motion.frequency,motion.sampling, motion.direction);
			}
			else if (motion.type=="randomwalk") {
				mwa.p_mw = new RandomMotionWalker(motion.displace,motion.seed,motion.sampling, motion.direction);				
			}
			else if (motion.type=="brownian") {
				mwa.p_mw = new BrownianMotionWalker(motion.displace,motion.seed,motion.sampling, motion.direction);				
			}
			else if (motion.type=="rotationalbrownian") {
				mwa.p_mw = new RotationalBrownianMotionWalker(motion.displace,motion.seed,motion.sampling);				
			}
			else if (motion.type=="localbrownian") {
				mwa.p_mw = new LocalBrownianMotionWalker(motion.radius,motion.displace,motion.seed,motion.sampling, motion.direction);							    
			}
			else if (motion.type=="none") {
                mwa.p_mw = NULL;
			}
			else {
                Err::Inst()->write("Motion type not understood");
                throw;
			}
			
			
			if (mwa.p_mw!=NULL) {
				m_motion_walkers.push_back(mwa);
			}
			
		}
	}
	
	// construct a helper for the alignments:
	if (Params::Inst()->sample.alignments.size()>0) {
		for(size_t i = 0; i < Params::Inst()->sample.alignments.size(); ++i)
		{
            CoordinateSetAlignment a;
            
            std::string type = Params::Inst()->sample.alignments[i].type;
            std::string order = Params::Inst()->sample.alignments[i].order;
            std::string selection = Params::Inst()->sample.alignments[i].selection;
            std::string rtype = Params::Inst()->sample.alignments[i].reference.type;
            std::string rselection = Params::Inst()->sample.alignments[i].reference.selection;    
            
            a.type = type;
            a.selection = selection;
            a.reference_selection = rselection;
            a.p_reference = NULL;

            // read and store reference now
            if (rtype=="frame") {
                size_t framenumber = Params::Inst()->sample.alignments[i].reference.frame;
                if (framenumber>=this->size()) {
                    Warn::Inst()->write("reference frame number in alignment larger than size of frameset. Setting to last frame!");
                    framenumber = this->size()-1;
                }
                Frame frame = frames.load(framenumber);
                a.p_reference = new CartesianCoordinateSet(frame,p_atoms->selections[a.reference_selection]);
            } else if (rtype=="file") {
                std::string file = Params::Inst()->sample.alignments[i].reference.file;
                std::string format = Params::Inst()->sample.alignments[i].reference.format;
                size_t framenumber = Params::Inst()->sample.alignments[i].reference.frame;
                if (format=="pdb") {
                    PDBFrameset fs(file,0);
                    fs.generate_index();
                    if (framenumber>=fs.number_of_frames) {
                        Warn::Inst()->write("reference frame number in alignment larger than size of frameset. Setting to last frame!");
                        framenumber = fs.number_of_frames-1;
                    }
                    Frame cf;		
                	fs.read_frame(framenumber,cf);
                	a.p_reference = new CartesianCoordinateSet(cf,p_atoms->selections[a.reference_selection]);
                } else {
                    Err::Inst()->write(string("File format for alignment reference not understood, format=") + format);
                    throw;
                }
	        } else if (rtype=="instant"){
				a.p_reference = NULL; // instructs code to use current coordinates
			} else {
                Err::Inst()->write(string("Reference type not understood:")+rtype );
                throw;
            }
            
            if (order=="pre") {
                m_prealignments.push_back(a);
            } else if (order=="post") {
                m_postalignments.push_back(a);                
            } else {
                Err::Inst()->write("Ordering of alignment not understood. Must be pre or post.");
                throw;
            }
		}
	}	
    
}

// std::vector<CartesianCoor3D> CoordinateSets::get_prealignmentvectors(size_t framenumber) {
//     return m_prealignmentvectors[framenumber];
// }
// 
// std::vector<CartesianCoor3D> CoordinateSets::get_postalignmentvectors(size_t framenumber) {
//     return m_postalignmentvectors[framenumber];    
// }

//void CoordinateSets::add_postalignment(std::string selection,std::string type) {
//    m_postalignments.push_back(make_pair(selection,type));
//}


void CoordinateSets::set_representation(CoordinateRepresentation representation) {
    m_representation = representation;
}

CoordinateRepresentation CoordinateSets::get_representation() {
    return m_representation;
}

void CoordinateSets::set_atoms(Atoms& atoms) {
    p_atoms = &atoms;
}

CoordinateSet* CoordinateSets::load(size_t framenumber) {
	
	Frame frame = frames.load(framenumber);
	CartesianCoordinateSet cset(frame,p_atoms->selections["system"]);

    if ( frame.x.size() != p_atoms->selections["system"]->size() ) {
        Err::Inst()->write(string("Wrong number of atoms, frame:")+boost::lexical_cast<string>(framenumber));
        Err::Inst()->write(string("Atoms in Frame: ")+boost::lexical_cast<string>(frame.x.size()));
        Err::Inst()->write(string("Atoms in Selection: ")+boost::lexical_cast<string>(p_atoms->selections["system"]->size()));
        throw;
    }

    // align here
    // m_prealignmentvectors[framenumber].clear();
    for(size_t i = 0; i < m_prealignments.size(); ++i)
    {
        CoordinateSetAlignment& a = m_prealignments[i];

		CoordinateSet* p_reference = NULL;
		IAtomselection* p_refselection=NULL;
		if (a.p_reference!=NULL) {
			p_reference=a.p_reference;
			p_refselection=p_atoms->selections[a.reference_selection];
		} else {
			p_reference=&cset;
			p_refselection=p_atoms->selections["system"];
		}

        if (a.type=="center") {
            CartesianCoor3D origin = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.reference_selection]);
            cset.translate(-1.0*origin,p_atoms->selections["system"], p_atoms->selections[a.selection]);
            // m_prealignmentvectors[framenumber].push_back(-1.0*origin);
        } else if (a.type=="fittrans") {
			CartesianCoor3D ref = CenterOfMass(*p_atoms,*p_reference,p_refselection);
            CartesianCoor3D pos = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection]);
            cset.translate(ref-pos,p_atoms->selections["system"], p_atoms->selections[a.selection]);
        } else if (a.type=="fitrottrans") {
            Fit(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection],*p_reference,p_refselection);            
        } else if (a.type=="fitrot") {
            CartesianCoor3D pos = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection]);
            Fit(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection],*p_reference,p_refselection);            
            // fit moves the center of mass, for rotational alignment only we have to correct the center of mass
            cset.translate(pos,p_atoms->selections["system"], p_atoms->selections[a.selection]);
        } else {
            Err::Inst()->write(string("Fitting routine not understood: ")+a.type);
            Err::Inst()->write(string("Use either of: center , fittrans , fitrottrans, fitrot"));
            throw;
        }
    }
        
	// apply motions, operate in cartesian space
	for(size_t i = 0; i < m_motion_walkers.size(); ++i)
	{
		MotionWalkerAlignment& mwa = m_motion_walkers[i];

		CartesianCoor3D ref(0,0,0);
		if (mwa.p_reference!=NULL) {
			ref = CenterOfMass(*p_atoms,*(mwa.p_reference),p_atoms->selections[mwa.reference_selection]);
		} else {
			ref = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[mwa.reference_selection]);
		}

//        CartesianCoor3D pos = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[mwa.selection]);
        cset.translate(-1.0*ref,p_atoms->selections["system"], p_atoms->selections[mwa.selection]);
		cset.transform(mwa.p_mw->transform(framenumber),p_atoms->selections["system"], p_atoms->selections[mwa.selection]);
        // fit moves the center of mass, for rotational alignment only we have to correct the center of mass
        cset.translate(ref,p_atoms->selections["system"], p_atoms->selections[mwa.selection]);
	}		

    // align here
    // m_postalignmentvectors[framenumber].clear();
    for(size_t i = 0; i < m_postalignments.size(); ++i)
    {
        CoordinateSetAlignment& a = m_postalignments[i];

		CoordinateSet* p_reference = NULL;
		IAtomselection* p_refselection=NULL;
		if (a.p_reference!=NULL) {
			p_reference=a.p_reference;
			p_refselection=p_atoms->selections[a.reference_selection];
		} else {
			p_reference=&cset;
			p_refselection=p_atoms->selections["system"];
		}

        if (a.type=="center") {
            CartesianCoor3D origin = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.reference_selection]);
            cset.translate(-1.0*origin,p_atoms->selections["system"], p_atoms->selections[a.selection]);
            // m_prealignmentvectors[framenumber].push_back(-1.0*origin);
        } else if (a.type=="fittrans") {
			CartesianCoor3D ref = CenterOfMass(*p_atoms,*p_reference,p_refselection);
            CartesianCoor3D pos = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection]);
            cset.translate(ref-pos,p_atoms->selections["system"], p_atoms->selections[a.selection]);
        } else if (a.type=="fitrottrans") {
            Fit(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection],*p_reference,p_refselection);            
        } else if (a.type=="fitrot") {
            CartesianCoor3D pos = CenterOfMass(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection]);
            Fit(*p_atoms,cset,p_atoms->selections["system"],p_atoms->selections[a.selection],*p_reference,p_refselection);            
            // fit moves the center of mass, for rotational alignment only we have to correct the center of mass
            cset.translate(pos,p_atoms->selections["system"], p_atoms->selections[a.selection]);
        } else {
            Err::Inst()->write(string("Fitting routine not understood: ")+a.type);
            Err::Inst()->write(string("Use either of: center , fittrans , fitrottrans, fitrot"));
            throw;
        }

    }    
    
    // reduce the coordinate set to the target selection
	CartesianCoordinateSet* pcset_reduced = new CartesianCoordinateSet(cset,p_atoms->selections["system"],p_selection);
    
    // convert to the current representation
	if (m_representation==CARTESIAN) {
        return pcset_reduced;
	} else if (m_representation==SPHERICAL) {
	    return (new SphericalCoordinateSet(*pcset_reduced) ); 	    
	} else if (m_representation==CYLINDRICAL) {
	    return (new CylindricalCoordinateSet(*pcset_reduced,Params::Inst()->scattering.average.orientation.axis) ); 	    
	} else {
        Err::Inst()->write("CoordinateSets::load: representation not understood");
        throw;
	}

}

void CoordinateSets::write_xyz(std::string filename) {
    ofstream ofile(filename.c_str());
	for(size_t i = 0; i < size(); ++i) {
        CoordinateSet* pcset = load(i);
        ofile << pcset->size() << endl;
        ofile << "generated by s_coordump" << endl;
        for(size_t j = 0; j < pcset->size(); ++j)
        {
            ofile << j << " " << pcset->c1[j] << " " << pcset->c2[j] << " " << pcset->c3[j] << endl;
        }
        delete pcset;
	}
    ofile.close();
}

void CoordinateSets::set_selection(IAtomselection* selection) {
	p_selection = selection;
}

IAtomselection& CoordinateSets::get_selection() {
	return *p_selection;
}


// end of file
