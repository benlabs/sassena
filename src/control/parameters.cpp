/** \file
This file contains the parameters class which contains the settings and values used to adjust the program control flow to the user input files.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "control/parameters.hpp"

// standard header
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// special library headers
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>	
#include <boost/thread.hpp>

// other headers
#include "exceptions/exceptions.hpp"
#include "log/log.hpp"
#include "io/xml_interface.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Params section
////////////////////////////////////////////////////////////////////////////////////////////////////


string Params::get_filepath(string fname) {
	using namespace boost::filesystem;

	path fpath(fname);
	string fdir;
	
	if (fpath.parent_path().is_complete()) 
		fdir = fpath.parent_path().string(); 
	else if (!config_rootpath.empty())
		fdir = (path(config_rootpath) / fpath.parent_path()).string();
	else 
		fdir = (initial_path() / fpath.parent_path()).string();
	return (path(fdir) / fpath.filename()).string();
}

void Params::write_xml_to_file(std::string filename) {
    std::ofstream conf(filename.c_str());
    
	conf << write_xml();

    conf.close();
}

void Params::read_xml(std::string filename) {
	
	// store the configuration in raw format into an internal buffer, 
	// this is for keeping history and allow for handing through data
	ifstream conffile(filename.c_str());
	while (conffile.good()) {
		char c=conffile.get();
		if (conffile.good()) {
			rawconfig.push_back(c);			
		}
	}
	conffile.close();
	
	// START OF sample section //	
	XMLInterface xmli(filename);
    
	Info::Inst()->write(string("Reading parameters from XML file: ")+filename);

	xmli.dump(config);
	// now read the parameters
	
	sample.structure.file="sample.pdb";
	sample.structure.filepath = get_filepath(sample.structure.file);
	sample.structure.format="pdb";
    
	if (xmli.exists("//sample")) {
    	if (xmli.exists("//sample/structure")) {
    	    if (xmli.exists("//sample/structure/file")) {
            	sample.structure.file   = xmli.get_value<std::string>("//sample/structure/file");
				sample.structure.filepath = get_filepath(sample.structure.file);
	        }
    	    if (xmli.exists("//sample/structure/format")) {
            	sample.structure.format   = xmli.get_value<std::string>("//sample/structure/format");
	        }
    	}   
    	if (xmli.exists("//sample/selections")) {
        	vector<XMLElement> selections = xmli.get("//sample/selections/selection");
        	for(size_t i = 0; i < selections.size(); ++i)
        	{
        		xmli.set_current(selections[i]);
                string type = "index";

                if (xmli.exists("./type")) type = xmli.get_value<string>("./type");
                
                if (type=="index") {
                    string sname = string("_")+boost::lexical_cast<string>(i); 
                    vector<size_t> ids;
                    
                    if (xmli.exists("./name")) sname = xmli.get_value<string>("./name") ;

                    vector<XMLElement> indexes = xmli.get("./index");
                    for(size_t pi = 0; pi < indexes.size(); ++pi)
                    {
                        xmli.set_current(indexes[pi]);
                        size_t index = xmli.get_value<size_t>(".");
                        ids.push_back(index);
                    }                        

                    sample.selections[sname] = new SampleIndexSelectionParameters(ids);         
                    Info::Inst()->write(string("Creating Index Atomselection: ")+sname+string(" , elements:")+boost::lexical_cast<string>(ids.size()));

                } else if (type=="range") {
                    string sname = string("_")+boost::lexical_cast<string>(i); 
                    size_t from = 0;
                    size_t to = 0;
                    
                    if (xmli.exists("./name")) sname = xmli.get_value<string>("./name") ;

                    if (xmli.exists("./from")) from = xmli.get_value<size_t>("./from") ;
                    if (xmli.exists("./to")) to = xmli.get_value<size_t>("./to") ;
                    
                    sample.selections[sname] = new SampleRangeSelectionParameters(from,to);
                    Info::Inst()->write(string("Creating Range Atomselection: ")+sname+string(" , from:")+boost::lexical_cast<string>(from) + string(", to: ")+boost::lexical_cast<string>(to));

                } else if (type=="lexical") {
                    string sname = string("_")+boost::lexical_cast<string>(i); 
                    string expression = "";
                    
                    if (xmli.exists("./name")) sname = xmli.get_value<string>("./name") ;

                    if (xmli.exists("./expression")) expression = xmli.get_value<string>("./expression") ;
                    
                    sample.selections[sname] = new SampleLexicalSelectionParameters(expression);
                    Info::Inst()->write(string("Creating Lexical Atomselection: ")+sname+string(" , expression:")+expression);                    
                    
                } else if (type=="file") {
                    string sname = string("_")+boost::lexical_cast<string>(i); // NOT used by ndx file format
                    string file = "selection.pdb";
                    string filepath = get_filepath(file);
                    string format = "pdb";
                    string selector = "beta";
                    string expression = "1|1\\.0|1\\.00";
                    
                    if (xmli.exists("./name")) sname = xmli.get_value<string>("./name") ;
                    
                    if (xmli.exists("./format")) format = xmli.get_value<string>("./format") ;

                    // this is a convenience overwrite for the ndx file format
                    if (format=="ndx") {
                        filename = "index.ndx";
                        selector = "name";
                        expression = ".*";
                    }
                    if (xmli.exists("./file")) {
						file = xmli.get_value<string>("./file") ;
						filepath = get_filepath(file) ;						
					}
                    
                    if (xmli.exists("./selector")) selector = xmli.get_value<string>("./selector") ;
                    if (xmli.exists("./expression")) expression = xmli.get_value<string>("./expression") ;
                    
                    sample.selections[sname] = new SampleFileSelectionParameters(filepath,format,selector,expression);                    
                    Info::Inst()->write(string("Creating File Atomselection: ")+sname+string(" , file: ")+filepath+string(", format: ")+format+string(", selector: ")+selector+string(", expression: ")+expression);
                }
        	}
	    }
	    if (xmli.exists("//sample/framesets")) {
        	size_t def_first=0;	
        	size_t def_last=0;	
        	bool def_last_set=false; 
        	size_t def_stride = 1;
            size_t def_clones = 1;
            string def_format = "dcd";
			string def_file="sample.dcd";
        	if (xmli.exists("//sample/framesets/first"))   def_first  = xmli.get_value<size_t>("//sample/framesets/first");
        	if (xmli.exists("//sample/framesets/last"))  { def_last   = xmli.get_value<size_t>("//sample/framesets/last"); def_last_set = true; }
        	if (xmli.exists("//sample/framesets/stride"))  def_stride = xmli.get_value<size_t>("//sample/framesets/stride");
        	if (xmli.exists("//sample/framesets/clones"))  def_clones = xmli.get_value<size_t>("//sample/framesets/clones");
        	if (xmli.exists("//sample/framesets/format"))   def_first  = xmli.get_value<size_t>("//sample/framesets/format");

        	// read in frame information

        	vector<XMLElement> framesets = xmli.get("//sample/framesets/frameset");
        	for(size_t i = 0; i < framesets.size(); ++i)
        	{
        		xmli.set_current(framesets[i]);
        		SampleFramesetParameters fset;	
        		fset.first = def_first;	fset.last = def_last; fset.last_set = def_last_set; fset.stride = def_stride;
                fset.format = def_format;
                fset.clones = def_clones;
				fset.file = def_file;
				fset.filepath = get_filepath(fset.file);
				
        		if (xmli.exists("./file")) {	
					fset.file = xmli.get_value<string>("./file");
					fset.filepath = get_filepath(fset.file);
				} 
                boost::filesystem::path index_path = fset.filepath;
				fset.index = (index_path.parent_path() / index_path.stem()).string() + string(".tnx");
				fset.index_default = true;
				
        		if (xmli.exists("./format"))  fset.format = xmli.get_value<string>("./format");
        		if (xmli.exists("./first"))   fset.first  = xmli.get_value<size_t>("./first");
        		if (xmli.exists("./last"))  { fset.last   = xmli.get_value<size_t>("./last"); fset.last_set = true; }
        		if (xmli.exists("./stride"))  fset.stride = xmli.get_value<size_t>("./stride");
        		if (xmli.exists("./clones"))  fset.clones = xmli.get_value<size_t>("./clones");
        		if (xmli.exists("./index"))  {
					fset.index = xmli.get_value<std::string>("./index");
					fset.indexpath = get_filepath(fset.index);
					fset.index_default = false;
				}
        		sample.framesets.push_back(fset);
        		Info::Inst()->write(string("Added frames from ")+fset.filepath+string(" using format: ")+fset.format);
        		Info::Inst()->write(string("Options: first=")+boost::lexical_cast<string>(fset.first)
        		    +string(", last=")+boost::lexical_cast<string>(fset.last)
        		    +string(", lastset=")+boost::lexical_cast<string>(fset.last_set)
        		    +string(", stride=")+boost::lexical_cast<string>(fset.stride));		
        	}
        }
        
	    if (xmli.exists("//sample/motions")) {
	    
	    	vector<XMLElement> motions = xmli.get("//sample/motions/motion");
	    	for(size_t i = 0; i < motions.size(); ++i)
	    	{
	    		xmli.set_current(motions[i]);
	    		SampleMotionParameters motion;	
	    		motion.type = "linear";	
	    		motion.displace = 0.0; 
	    		motion.direction=CartesianCoor3D(1,0,0);
	    		motion.selection = "system";
	    		motion.seed = 0;
	    		motion.sampling = 1;			
	    		motion.frequency=0.001; // corresponds to one full cycle per 1000 frames, used for linear oscillation and rotation
	
	            motion.reference.type = "instant";
                motion.reference.frame = 0;
                motion.reference.file = sample.structure.file;
                motion.reference.filepath = get_filepath(motion.reference.file);
                motion.reference.format = sample.structure.format;
                motion.reference.selection = motion.selection;
    
	    		if (xmli.exists("./type"))   motion.type  = xmli.get_value<string>("./type");
	    		if (xmli.exists("./displace"))  motion.displace   = xmli.get_value<double>("./displace");
	    		if (xmli.exists("./frequency"))  motion.frequency   = xmli.get_value<double>("./frequency");			
	    		if (xmli.exists("./seed"))  motion.seed   = xmli.get_value<unsigned long>("./seed");			
	    		if (xmli.exists("./sampling"))  motion.sampling   = xmli.get_value<long>("./sampling");			
	    		if (xmli.exists("./selection"))  motion.selection   = xmli.get_value<string>("./selection");			
	    		if (xmli.exists("./direction")) {
	    			if (xmli.exists("./direction/x")) motion.direction.x   = xmli.get_value<double>("./direction/x");
	    			if (xmli.exists("./direction/y")) motion.direction.y   = xmli.get_value<double>("./direction/y");
	    			if (xmli.exists("./direction/z")) motion.direction.z   = xmli.get_value<double>("./direction/z");				
	    		} 
                motion.radius = motion.displace*10;
	    		if (xmli.exists("./radius"))  motion.radius   = xmli.get_value<double>("./radius");			
        
				
    	    	if (xmli.exists("./reference")) {
        	    	if (xmli.exists("./reference/type")) motion.reference.type   = xmli.get_value<string>("./reference/type");
        	    	if (xmli.exists("./reference/frame")) motion.reference.frame   = xmli.get_value<size_t>("./reference/frame");
        	    	if (xmli.exists("./reference/file")) {
						motion.reference.file   = xmli.get_value<string>("./reference/file");
						motion.reference.filepath   = get_filepath(motion.reference.file);						
					}
        	    	if (xmli.exists("./reference/format")) motion.reference.format   = xmli.get_value<string>("./reference/format");
        	    	if (xmli.exists("./reference/selection")) motion.reference.selection   = xmli.get_value<string>("./reference/selection");
	    	    }
	    	    if (motion.type=="file") {
        	    	Info::Inst()->write(string("Reference for motion alignment: type=")+motion.reference.type+string(", file=")+motion.reference.filepath+string(", format=")+motion.reference.format+string(", frame=")+boost::lexical_cast<string>(motion.reference.frame));	    	        
	    	    } else if (motion.type=="frame"){
        	    	Info::Inst()->write(string("Reference for motion alignment: type=")+motion.reference.type+string(", frame=")+boost::lexical_cast<string>(motion.reference.frame));
        	    	Info::Inst()->write(string("Motion Alignment uses unprocessed coordinates (No alignment and no applied motions)"));
	    	    } else if (motion.type=="instant") {
        	    	Info::Inst()->write(string("Instant motion alignment (Uses coordinates within the current frame)"));		
				}
	    	    Info::Inst()->write(string("Selection used for alignment with reference=")+motion.reference.selection);	    	            	    	

	    		sample.motions.push_back(motion);
	    		Info::Inst()->write(string("Adding additional motion to sample: type=")+motion.type
	    		+string(", displacement=")+boost::lexical_cast<string>(motion.displace)
	    		+string(", selection=")+motion.selection);
	    	}
	    }
	    
        if (xmli.exists("//sample/alignments")) {

    	    vector<XMLElement> alignments = xmli.get("//sample/alignments/alignment");
    	    for(size_t i = 0; i < alignments.size(); ++i)
    	    {
    	    	xmli.set_current(alignments[i]);
    	    	SampleAlignmentParameters alignment;	
    	    	alignment.type = "center";	 
    	    	alignment.selection = "system";
                alignment.order = "pre";
    	    	if (xmli.exists("./type"))   alignment.type  = xmli.get_value<string>("./type");
    	    	if (xmli.exists("./selection")) alignment.selection   = xmli.get_value<string>("./selection");
    	    	if (xmli.exists("./order"))  alignment.order   = xmli.get_value<string>("./order");			

    	    	Info::Inst()->write(string("Adding additional alignment to sample: type=")+alignment.type+string(", selection=")+alignment.selection+string(", order=")+alignment.order);

                alignment.reference.type = "frame";
                alignment.reference.frame = 0;
                alignment.reference.file = sample.structure.file;
				alignment.reference.filepath   = get_filepath(alignment.reference.file);
                alignment.reference.format = sample.structure.format;
                alignment.reference.selection = alignment.selection;
                
    	    	if (xmli.exists("./reference")) {
        	    	if (xmli.exists("./reference/type")) alignment.reference.type   = xmli.get_value<string>("./reference/type");
        	    	if (xmli.exists("./reference/frame")) alignment.reference.frame   = xmli.get_value<size_t>("./reference/frame");
        	    	if (xmli.exists("./reference/file")) {
						alignment.reference.file   = xmli.get_value<string>("./reference/file");
						alignment.reference.filepath   = get_filepath(alignment.reference.file);						
					}
        	    	if (xmli.exists("./reference/format")) alignment.reference.format   = xmli.get_value<string>("./reference/format");
        	    	if (xmli.exists("./reference/selection")) alignment.reference.selection   = xmli.get_value<string>("./reference/selection");
	    	    }
	    	    if (alignment.type=="file") {
        	    	Info::Inst()->write(string("Reference for alignment: type=")+alignment.reference.type+string(", file=")+alignment.reference.filepath+string(", format=")+alignment.reference.format+string(", frame=")+boost::lexical_cast<string>(alignment.reference.frame));	    	        
	    	    } else if (alignment.type=="frame"){
        	    	Info::Inst()->write(string("Reference for alignment: type=")+alignment.reference.type+string(", frame=")+boost::lexical_cast<string>(alignment.reference.frame));
        	    	Info::Inst()->write(string("Alignment uses unprocessed coordinates (No alignment and no applied motions)"));
	    	    }
	    	    Info::Inst()->write(string("Selection used for alignment with reference=")+alignment.reference.selection);	    	            	    	

    	    	sample.alignments.push_back(alignment);
    	    }
        }	
	
	}	
	
    stager.dump=false;
    stager.file="dump.dcd";
    stager.filepath=get_filepath(stager.file);
    stager.format="dcd";
    stager.target = "system";
    stager.mode = "frames";

    if (xmli.exists("//stager")) {
    	if (xmli.exists("//stager/dump")) {
    		stager.dump = xmli.get_value<bool>("//stager/dump");
    	}
    	if (xmli.exists("//stager/file")) {
    		stager.file = xmli.get_value<string>("//stager/file");
			stager.filepath=get_filepath(stager.file);
    	}
    	if (xmli.exists("//stager/format")) {
    		stager.format = xmli.get_value<string>("//stager/format");
    	}
    	if (xmli.exists("//stager/target")) {
    		stager.target = xmli.get_value<string>("//stager/target");        	
    	}
    	if (xmli.exists("//stager/mode")) {
    		stager.mode = xmli.get_value<string>("//stager/mode");
    	}
    }
	Info::Inst()->write(string("stager.target=")+stager.target);		
	
	scattering.type="all";
	scattering.background.factor = 0.0;

    scattering.dsp.type="autocorrelate";
    scattering.dsp.method="fftw";
	// defaults
	scattering.average.orientation.type = "none";
	scattering.average.orientation.axis = CartesianCoor3D(0,0,1);
	scattering.average.orientation.vectors.type = "sphere";
	scattering.average.orientation.vectors.algorithm = "boost_uniform_on_sphere";
	scattering.average.orientation.vectors.resolution = 100;
	scattering.average.orientation.vectors.seed = 0;
	scattering.average.orientation.vectors.file = "qvector-orientations.txt";
	scattering.average.orientation.vectors.filepath = get_filepath(scattering.average.orientation.vectors.file);
	scattering.average.orientation.multipole.type = "sphere";
	scattering.average.orientation.multipole.moments.resolution = 20;
    scattering.average.orientation.multipole.moments.file = "moments.txt";
    scattering.average.orientation.multipole.moments.filepath = get_filepath(scattering.average.orientation.multipole.moments.file);
	scattering.average.orientation.exact.type = "sphere";

    scattering.signal.file = "signal.h5";
    scattering.signal.filepath = get_filepath(scattering.signal.file);
    scattering.signal.fqt = true;
    scattering.signal.fq0 = true;
    scattering.signal.fq = true;
    scattering.signal.fq2 = true;

    if (xmli.exists("//scattering")) {
    	if (xmli.exists("//scattering/type")) {
    		scattering.type = xmli.get_value<string>("//scattering/type");
    	}

    	if (xmli.exists("//scattering/target")) {
			Err::Inst()->write("scattering.target is obsolete. Use stager.target instead.");
			throw;
    	}

    	Info::Inst()->write(string("scattering.type=")+scattering.type);

    	if (xmli.exists("//scattering/background")) {
    		if (xmli.exists("//scattering/background/factor")) {
    		    scattering.background.factor = xmli.get_value<double>("//scattering/background/factor");
    		    Info::Inst()->write(string("scattering.background.factor=")+boost::lexical_cast<string>(scattering.background.factor));		    		    
    	    }

            if (xmli.exists("//scattering/background/kappas")) {

        	    vector<XMLElement> kappas = xmli.get("//scattering/background/kappas/kappa");
        	    for(size_t i = 0; i < kappas.size(); ++i)
        	    {
        	    	xmli.set_current(kappas[i]);
        	    	ScatteringBackgroundKappaParameters kappa;	
        	    	kappa.selection = "system";	 
        	    	kappa.value = 1.0;
        	    	if (xmli.exists("./selection"))   kappa.selection  = xmli.get_value<string>("./selection");
        	    	if (xmli.exists("./value"))  kappa.value   = xmli.get_value<double>("./value");			

        	    	scattering.background.kappas.push_back(kappa);
        	    	Info::Inst()->write(string("Rescaling volumes: selection=")+kappa.selection+string(", value=")+boost::lexical_cast<string>(kappa.value));
        	    }
            }
    	}	
    	// generating qqqvectors, i.e. the spectrum
	    if (xmli.exists("//scattering/vectors")) {
            
            string vt = "single";
            
            if (xmli.exists("//scattering/vectors/type"))  vt = xmli.get_value<string>("//scattering/vectors/type");
    	    if (vt=="single") {	
    	        double x = 1.0;
                double y = 0;
                double z = 0;
                
                if (xmli.exists("//scattering/vectors/single/x")) {
        	    	if (xmli.exists("//scattering/vectors/single/x")) x = xmli.get_value<double>("//scattering/vectors/single/x");
        	    	if (xmli.exists("//scattering/vectors/single/y")) y = xmli.get_value<double>("//scattering/vectors/single/y");
        	    	if (xmli.exists("//scattering/vectors/single/z")) z = xmli.get_value<double>("//scattering/vectors/single/z");                    
                }
            
    	    	scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
    	    }
    	    else if (vt=="scans") {	
                
                if (xmli.exists("//scattering/vectors/scans")) {
    	    	    vector<XMLElement> scans = xmli.get("//scattering/vectors/scans/scan");
            
    	    	    for(size_t i = 0; i < scans.size(); ++i)
    	    	    {
    	    		    xmli.set_current(scans[i]);
    	    		    ScatteringVectorsScanParameters sc;
            	        sc.basevector.x = 1.0;
                        sc.basevector.y = 0;
                        sc.basevector.z = 0;
                        sc.from = 0;
                        sc.to = 1;
                        sc.points = 100;
                        sc.exponent = 1.0;
                        
                        if (xmli.exists("./base")) {
        	    		    if (xmli.exists("./base/x"))         sc.basevector.x = xmli.get_value<double>("./base/x");
        	    		    if (xmli.exists("./base/y"))         sc.basevector.y = xmli.get_value<double>("./base/y");
        	    		    if (xmli.exists("./base/z"))         sc.basevector.z = xmli.get_value<double>("./base/z");                            
                        }
    	    		    if (xmli.exists("./from"))      sc.from     = xmli.get_value<double>("./from");
    	    		    if (xmli.exists("./to"))        sc.to       = xmli.get_value<double>("./to");
    	    		    if (xmli.exists("./points"))    sc.points   = xmli.get_value<size_t>("./points");
    	    		    if (xmli.exists("./exponent"))  sc.exponent = xmli.get_value<double>("./exponent");
    	    		    scattering.qvectors.scans.push_back(sc);
    	    	    }
	    	    }
            
    	    	scattering.qvectors.create_from_scans();
    	    }
    	    else if (vt=="file") {
                string qqqfilename = "qvectors.txt";
                if (xmli.exists("//scattering/vectors/file")) qqqfilename = get_filepath(xmli.get_value<string>("//scattering/vectors/file"));
    	    	ifstream qqqfile(qqqfilename.c_str());
            
    	    	double x,y,z; 
    	    	while (qqqfile >> x >> y >> z) {
    	    		scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
    	    	}
    	    }
        }
        
        if (scattering.qvectors.size()==0) {
            Err::Inst()->write("No q vectors generated. Check the scattering.vectors section for errors.");
            throw;
        }

    	if (xmli.exists("//scattering/dsp")) {
    		if (xmli.exists("//scattering/dsp/type")) {
    			scattering.dsp.type = xmli.get_value<string>("//scattering/dsp/type");
    		    Info::Inst()->write(string("scattering.dsp.type=")+scattering.dsp.type);
    		}
    		if (xmli.exists("//scattering/dsp/method")) {
    			scattering.dsp.method = xmli.get_value<string>("//scattering/dsp/method");
    		    Info::Inst()->write(string("scattering.dsp.method=")+scattering.dsp.method);
    		}
    	}
    	
	    if (xmli.exists("//scattering/average")) {
	    	if (xmli.exists("//scattering/average/orientation")) {
	    	    if (xmli.exists("//scattering/average/orientation/axis")) { // count vectors ... , or order for multipole...
    				scattering.average.orientation.axis.x = xmli.get_value<double>("//scattering/average/orientation/axis/x");
    				scattering.average.orientation.axis.y = xmli.get_value<double>("//scattering/average/orientation/axis/y");
    				scattering.average.orientation.axis.z = xmli.get_value<double>("//scattering/average/orientation/axis/z");					
    			}									
	    	    
	    		if (xmli.exists("//scattering/average/orientation/type")) { // vectors multipole
	    			scattering.average.orientation.type = xmli.get_value<string>("//scattering/average/orientation/type");
	    		}
	    	    Info::Inst()->write(string("scattering.average.orientation.type=")+scattering.average.orientation.type);
	    		
	    		if (scattering.average.orientation.type=="vectors") {
	    			if (xmli.exists("//scattering/average/orientation/vectors/type")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.type = xmli.get_value<string>("//scattering/average/orientation/vectors/type");
                        Info::Inst()->write(string("scattering.average.orientation.vectors.type=")+scattering.average.orientation.vectors.type);				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/vectors/algorithm")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.algorithm = xmli.get_value<string>("//scattering/average/orientation/vectors/algorithm");
                        Info::Inst()->write(string("scattering.average.orientation.vectors.algorithm=")+scattering.average.orientation.vectors.algorithm);				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/vectors/resolution")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.resolution = xmli.get_value<long>("//scattering/average/orientation/vectors/resolution");
                        Info::Inst()->write(string("scattering.average.orientation.vectors.resolution=")+boost::lexical_cast<string>(scattering.average.orientation.vectors.resolution));				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/vectors/seed")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.seed = xmli.get_value<unsigned long>("//scattering/average/orientation/vectors/seed");
                        Info::Inst()->write(string("scattering.average.orientation.vectors.seed=")+boost::lexical_cast<string>(scattering.average.orientation.vectors.seed));
	    			}	
	    			if (xmli.exists("//scattering/average/orientation/vectors/file")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.file = xmli.get_value<string>("//scattering/average/orientation/vectors/file");
	    				scattering.average.orientation.vectors.filepath = get_filepath(scattering.average.orientation.vectors.file);
	                    Info::Inst()->write(string("scattering.average.orientation.vectors.file=")+boost::lexical_cast<string>(scattering.average.orientation.vectors.file));				    
                        Info::Inst()->write(string("scattering.average.orientation.vectors.filepath=")+boost::lexical_cast<string>(scattering.average.orientation.vectors.filepath));				    
	    			}		
	    				    			
	    			scattering.average.orientation.vectors.create();
	    			
	    		} else if (scattering.average.orientation.type=="multipole") {			    
	    			if (xmli.exists("//scattering/average/orientation/multipole/type")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.multipole.type = xmli.get_value<string>("//scattering/average/orientation/multipole/type");
                        Info::Inst()->write(string("scattering.average.orientation.multipole.type=")+scattering.average.orientation.multipole.type);				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/multipole/moments")) { // count vectors ... , or order for multipole...
    	    			if (xmli.exists("//scattering/average/orientation/multipole/moments/type")) { // count vectors ... , or order for multipole...
	    				    scattering.average.orientation.multipole.moments.type = xmli.get_value<string>("//scattering/average/orientation/multipole/moments/type");
                            Info::Inst()->write(string("scattering.average.orientation.multipole.moments.type=")+boost::lexical_cast<string>(scattering.average.orientation.multipole.moments.type));
                        }
    	    			if (xmli.exists("//scattering/average/orientation/multipole/moments/resolution")) { // count vectors ... , or order for multipole...
	    				    scattering.average.orientation.multipole.moments.resolution = xmli.get_value<long>("//scattering/average/orientation/multipole/moments/resolution");
                            Info::Inst()->write(string("scattering.average.orientation.multipole.moments.resolution=")+boost::lexical_cast<string>(scattering.average.orientation.multipole.moments.resolution));				    
                        }
    	    			if (xmli.exists("//scattering/average/orientation/multipole/moments/file")) { // count vectors ... , or order for multipole...
	    				    scattering.average.orientation.multipole.moments.file = xmli.get_value<string>("//scattering/average/orientation/multipole/moments/file");
	    				    scattering.average.orientation.multipole.moments.filepath = get_filepath(scattering.average.orientation.multipole.moments.file);
                            Info::Inst()->write(string("scattering.average.orientation.multipole.moments.file=")+boost::lexical_cast<string>(scattering.average.orientation.multipole.moments.file));				    
                            Info::Inst()->write(string("scattering.average.orientation.multipole.moments.filepath=")+boost::lexical_cast<string>(scattering.average.orientation.multipole.moments.filepath));				    
                        }
	    			}
	    			
	    			scattering.average.orientation.multipole.moments.create();
	    		} else if (scattering.average.orientation.type!="none") {
	    			Err::Inst()->write(string("Orientation averaging type not understood: ")+scattering.average.orientation.type);
	    			throw;
	    		}
	    	}
	    }
	    
	    if (xmli.exists("//scattering/signal")) {
	        if (xmli.exists("//scattering/signal/file")) {
        		scattering.signal.file = xmli.get_value<string>("//scattering/signal/file");	        
        		scattering.signal.filepath = get_filepath(scattering.signal.file);	        
	        } 
	        if (xmli.exists("//scattering/signal/fqt")) {
        		scattering.signal.fqt = xmli.get_value<bool>("//scattering/signal/fqt");	        
	        } 
	        if (xmli.exists("//scattering/signal/fq0")) {
        		scattering.signal.fq0 = xmli.get_value<bool>("//scattering/signal/fq0");	        
	        }
	        if (xmli.exists("//scattering/signal/fq")) {
        		scattering.signal.fq = xmli.get_value<bool>("//scattering/signal/fq");	        
	        }
	        if (xmli.exists("//scattering/signal/fq2")) {
        		scattering.signal.fq2 = xmli.get_value<bool>("//scattering/signal/fq2");	        
	        }	         
	    }
	    Info::Inst()->write(string("scattering.signal.file=")+scattering.signal.file);	
	    Info::Inst()->write(string("scattering.signal.filepath=")+scattering.signal.filepath);	
	    Info::Inst()->write(string("scattering.signal.fqt=")+boost::lexical_cast<string>(scattering.signal.fqt));	
	    Info::Inst()->write(string("scattering.signal.fq0=")+boost::lexical_cast<string>(scattering.signal.fq0));	            
	    Info::Inst()->write(string("scattering.signal.fq=")+boost::lexical_cast<string>(scattering.signal.fq));	
	    Info::Inst()->write(string("scattering.signal.fq2=")+boost::lexical_cast<string>(scattering.signal.fq2));	    
    }  
	// END OF scattering section //

	// START OF limits section //

    // assign default memory limits:    
    limits.stage.memory.buffer = 100*1024*1024;
    limits.stage.memory.data   = 500*1024*1024;
    
    limits.signal.chunksize = 10000;

    limits.computation.cores = 1;
    limits.computation.processes = 1;
    limits.computation.threads = 1;
    limits.computation.memory.result_buffer = 100*1024*1024; // 100MB 
    limits.computation.memory.signal_buffer = 100*1024*1024; // 100MB 
    limits.computation.memory.exchange_buffer = 100*1024*1024; // 100MB 
    limits.computation.memory.alignpad_buffer = 200*1024*1024; // 200MB 
    limits.computation.memory.scale = 1;

    limits.services.signal.memory.server = 100*1024*1024; // 100MB
    limits.services.signal.memory.client = 10*1024*1024; // 10MB
    limits.services.signal.times.serverflush = 600; // 600 seconds
    limits.services.signal.times.clientflush = 600; // 600 seconds

    limits.services.monitor.delay = 1; // 1 second
    limits.services.monitor.sampling = 0; // 0 = automatic 
    
    limits.decomposition.utilization = 0.95; // 5% max loss
    limits.decomposition.partitions.automatic = true; // pick number of independent partitions based on some heuristics
    limits.decomposition.partitions.size = 1; // not used if automatic = true, if false -> this determines the partition size

	if (xmli.exists("//limits")) {       	
	    if (xmli.exists("//limits/stage")) {	         
        	if (xmli.exists("//limits/stage/memory")) {
            	if (xmli.exists("//limits/stage/memory/buffer")) {
        	        limits.stage.memory.buffer = xmli.get_value<size_t>("//limits/stage/memory/buffer");
    	        }	        
            	if (xmli.exists("//limits/stage/memory/data")) {
        	        limits.stage.memory.data = xmli.get_value<size_t>("//limits/stage/memory/data");
    	        }	            	        
	        }	        
	    }
	    
	    if (xmli.exists("//limits/signal")) {
    	    if (xmli.exists("//limits/signal/chunksize")) {
    	        limits.signal.chunksize = xmli.get_value<size_t>("//limits/signal/chunksize");
            }                  
	    }

    	if (xmli.exists("//limits/computation")) {
            size_t cores = boost::thread::hardware_concurrency();
        	if (xmli.exists("//limits/computation/cores")) {
    	        limits.computation.cores = xmli.get_value<size_t>("//limits/computation/cores");                
    	        Info::Inst()->write(string("Assuming ")+boost::lexical_cast<string>(limits.computation.cores)+string(" cores per machine"));
            } else {                
                Info::Inst()->write(string("Detect: Number of Processors per machine: ")+boost::lexical_cast<string>(cores));   
                limits.computation.cores = cores;
            }
        	if (xmli.exists("//limits/computation/processes")) {
    	        limits.computation.processes = xmli.get_value<size_t>("//limits/computation/processes");
	        }

        	if (xmli.exists("//limits/computation/threads")) {
    	        limits.computation.threads = xmli.get_value<size_t>("//limits/computation/threads");
	        }
        	if (xmli.exists("//limits/computation/memory")) {
            	if (xmli.exists("//limits/computation/memory/result_buffer")) {
    	            limits.computation.memory.result_buffer = xmli.get_value<size_t>("//limits/computation/memory/result_buffer");
    	        }
            	if (xmli.exists("//limits/computation/memory/signal_buffer")) {
    	            limits.computation.memory.signal_buffer = xmli.get_value<size_t>("//limits/computation/memory/signal_buffer");
    	        }
            	if (xmli.exists("//limits/computation/memory/exchange_buffer")) {
    	            limits.computation.memory.exchange_buffer = xmli.get_value<size_t>("//limits/computation/memory/exchange_buffer");
    	        }
            	if (xmli.exists("//limits/computation/memory/alignpad_buffer")) {
    	            limits.computation.memory.alignpad_buffer = xmli.get_value<size_t>("//limits/computation/memory/alignpad_buffer");
    	        }
            	if (xmli.exists("//limits/computation/memory/scale")) {
    	            limits.computation.memory.scale = xmli.get_value<size_t>("//limits/computation/memory/scale");
    	        }

	        }
	    }
	    
	    if (xmli.exists("//limits/services")) {
        	if (xmli.exists("//limits/services/signal")) {
            	if (xmli.exists("//limits/services/signal/memory")) {
                	if (xmli.exists("//limits/services/signal/memory/server")) {
        	            limits.services.signal.memory.server = xmli.get_value<size_t>("//limits/services/signal/memory/server");
    	            }
                	if (xmli.exists("//limits/services/signal/memory/client")) {
        	            limits.services.signal.memory.client = xmli.get_value<size_t>("//limits/services/signal/memory/client");
    	            }
    	        }
            	if (xmli.exists("//limits/services/signal/times")) {
                	if (xmli.exists("//limits/services/signal/times/serverflush")) {
        	            limits.services.signal.times.serverflush = xmli.get_value<size_t>("//limits/services/signal/times/serverflush");
    	            }
                	if (xmli.exists("//limits/services/signal/times/clientflush")) {
        	            limits.services.signal.times.clientflush = xmli.get_value<size_t>("//limits/services/signal/times/clientflush");
    	            }
    	        }    	        
	        }
	    }
	    if (xmli.exists("//limits/services")) {
        	if (xmli.exists("//limits/services/monitor")) {
            	if (xmli.exists("//limits/services/monitor/delay")) {
    	            limits.services.monitor.delay = xmli.get_value<size_t>("//limits/services/monitor/delay");
    	        }    	        
            	if (xmli.exists("//limits/services/monitor/sampling")) {
    	            limits.services.monitor.sampling = xmli.get_value<size_t>("//limits/services/monitor/sampling");
    	        }    	        
	        }
	    }	    
	    
    	if (xmli.exists("//limits/decomposition")) {
        	if (xmli.exists("//limits/decomposition/utilization")) {
			    limits.decomposition.utilization = xmli.get_value<double>("//limits/decomposition/utilization");
            }    	    
        	if (xmli.exists("//limits/decomposition/partitions")) {
            	if (xmli.exists("//limits/decomposition/partitions/automatic")) {
    			    limits.decomposition.partitions.automatic = xmli.get_value<bool>("//limits/decomposition/partitions/automatic");
                }
            	if (xmli.exists("//limits/decomposition/partitions/size")) {
    			    limits.decomposition.partitions.size = xmli.get_value<size_t>("//limits/decomposition/partitions/size");                    
                }
            }
        }		
	}

	// END OF limits section //
	// START OF debug section //

	debug.timer = false; // this adds a log message when a timer is started/stopped
	debug.barriers = false; // this de-/activates collective barriers before each collective operation, this way all nodes are synchronized before the communication takes place. This is an important step towards analysis of timing.
    debug.monitor.update = true;
    debug.iowrite.write = true;
    debug.iowrite.buffer = true;
    debug.print.orientations = false;
	if (xmli.exists("//debug")) {
		if (xmli.exists("//debug/timer")) {
			debug.timer = xmli.get_value<bool>("//debug/timer");
		}
		if (xmli.exists("//debug/barriers")) {
			debug.barriers = xmli.get_value<bool>("//debug/barriers");
		}
		if (xmli.exists("//debug/monitor")) {
    		if (xmli.exists("//debug/monitor/update")) {
		    	debug.monitor.update = xmli.get_value<bool>("//debug/monitor/update");
	    	}
		}
		if (xmli.exists("//debug/iowrite")) {
    		if (xmli.exists("//debug/iowrite/write")) {
		    	debug.iowrite.write = xmli.get_value<bool>("//debug/iowrite/write");
		    }
    		if (xmli.exists("//debug/iowrite/buffer")) {
		    	debug.iowrite.buffer = xmli.get_value<bool>("//debug/iowrite/buffer");
		    }
		}
		
		if (xmli.exists("//debug/print")) {
		    if (xmli.exists("//debug/print/orientations")) {
		    	debug.print.orientations = xmli.get_value<bool>("//debug/print/orientations");
	    	}
		}
		
	}
	
    database.type = "file";
    database.file = "db.xml";
    database.filepath = get_filepath(database.file);
    database.format = "xml";
	if (xmli.exists("//database")) {
	    if (xmli.exists("//database/type")) {
	        database.type = xmli.get_value<string>("//database/type");
        }
	    if (xmli.exists("//database/file")) {
	        database.file = xmli.get_value<string>("//database/file");
	        database.filepath = get_filepath(database.file);
        }
	    if (xmli.exists("//database/format")) {
	        database.format = xmli.get_value<string>("//database/format");
        }        
    }
};


boost::program_options::options_description Params::options() {

    namespace po = boost::program_options;
    po::options_description all("Overwrite options");
    
    po::options_description generic("Generic options");
    generic.add_options()
        ("help", "produce this help message")
        ("config", po::value<string>()->default_value("scatter.xml"),  "name of the xml configuration file")
    ;

    po::options_description sample("Sample related options");
    sample.add_options()
        ("sample.structure.file",po::value<string>()->default_value("sample.pdb"), "Structure file name")
        ("sample.structure.format",po::value<string>()->default_value("pdb"), "Structure file format")
    ;

    po::options_description stager("Data staging related options");
    stager.add_options()
        ("stager.target",po::value<string>()->default_value("system"), "Atom selection producing the signal (must be defined)")
        ("stager.dump",po::value<bool>()->default_value(false), "Do/Don't dump the postprocessed coordinates to a file")        
        ("stager.file",po::value<string>()->default_value("dump.dcd"), "Name of dump file")        
        ("stager.format",po::value<string>()->default_value("dcd"), "Format of dump file")        
    ;


    po::options_description scattering("Scattering related options");
    scattering.add_options()
        ("scattering.signal.file",po::value<string>()->default_value("signal.h5"), "name of the signal file")
    ;

    po::options_description limits("Computational limits related options");
    limits.add_options()
        ("limits.computation.threads",po::value<int>()->default_value(1), "Number of worker threads per process")
    ;
    
    all.add(generic);
    all.add(sample);
    all.add(stager);
    all.add(scattering);
    all.add(limits);

    return all;
}

void Params::overwrite_options(boost::program_options::variables_map& vm) {
    if (!vm["sample.structure.file"].defaulted()) {
        std::string val = vm["sample.structure.file"].as<string>();
        Info::Inst()->write(string("OVERWRITE sample.structure.file=")+val);
        Params::Inst()->sample.structure.file=val;
    }
    if (!vm["sample.structure.format"].defaulted()) {
        std::string val = vm["sample.structure.format"].as<string>();
        Info::Inst()->write(string("OVERWRITE sample.structure.format=")+val);
        Params::Inst()->sample.structure.format=val;
    }
    if (!vm["stager.target"].defaulted()) {
        std::string val = vm["stager.target"].as<string>();
        Info::Inst()->write(string("OVERWRITE stager.target=")+val);
        Params::Inst()->stager.target=val;
    }
    if (!vm["stager.dump"].defaulted()) {
        bool val = vm["stager.dump"].as<bool>();
        Info::Inst()->write(string("OVERWRITE stager.dump=")+boost::lexical_cast<string>(val));
        Params::Inst()->stager.dump=val;
    }
    if (!vm["stager.file"].defaulted()) {
        std::string val = vm["stager.file"].as<string>();
        Info::Inst()->write(string("OVERWRITE stager.file=")+val);
        Params::Inst()->stager.file=val;
    }
    if (!vm["stager.format"].defaulted()) {
        std::string val = vm["stager.format"].as<string>();
        Info::Inst()->write(string("OVERWRITE stager.format=")+val);
        Params::Inst()->stager.format=val;
    }

    if (!vm["scattering.signal.file"].defaulted()) {
        std::string val = vm["scattering.signal.file"].as<string>();
        Info::Inst()->write(string("OVERWRITE scattering.signal.file=")+val);
        Params::Inst()->scattering.signal.file=val;
    }
    if (!vm["limits.computation.threads"].defaulted()) {
        int val = vm["limits.computation.threads"].as<int>();
        Info::Inst()->write(string("OVERWRITE limits.computation.threads=")+boost::lexical_cast<string>(val));
        Params::Inst()->limits.computation.threads=val;
    }
}

void Params::init(int argc,char** argv) {

    namespace po = boost::program_options;
    po::variables_map vm;
    po::options_description ops = options();
    
    po::store(po::parse_command_line(argc, argv, ops), vm);
    po::notify(vm);
    
    if (vm.find("help")!=vm.end()) {
        cout << ops << endl;
        throw sassena::terminate_request();
    }

    if (vm["config"].defaulted()) {
        Info::Inst()->write("No configuration file specified. Will use default");
    }

    std::string filename = vm["config"].as<string>();
    if (boost::filesystem::exists(boost::filesystem::path(filename))) {
	    // make the directory of the main configuration file the root for all others
	    if (boost::filesystem::path(filename).is_complete()) 
		    config_rootpath = boost::filesystem::path(filename).parent_path().string();
	    else 
		    config_rootpath = ( boost::filesystem::initial_path() / boost::filesystem::path(filename).parent_path() ).string();		

	    Info::Inst()->write(string("Looking for configuration file: ") + filename);
	
	    read_xml(filename);
    } else {
        Err::Inst()->write(vm["config"].as<string>()+string(" does not exist!"));            
        throw;
    }
    
	
	Info::Inst()->write(string("Analyzing command line for parameter overwrites"));
	
    overwrite_options(vm);
}

SampleParameters::~SampleParameters() {
    for (std::map<std::string,SampleSelectionParameters*>::iterator i=selections.begin();i!=selections.end();i++) {
        delete i->second;
    }
}

void ScatteringAverageOrientationVectorsParameters::create() {
	
	if (type=="file") {
		Info::Inst()->write(string("Reading orientations for orientational averaging from file: ")+filepath);
        
		ifstream qqqfile(filepath.c_str());
	
		double x,y,z; 
		while (qqqfile >> x >> y >> z) {
			CartesianCoor3D q(x,y,z); 
			double ql = q.length();
			if (ql!=0) q = (1.0/ql) * q;
			this->push_back(q);
		}
		
	} else if (type=="sphere") {
		if (algorithm=="boost_uniform_on_sphere") {
	
			boost::mt19937 rng; // that's my random number generator
			rng.seed(seed);
			boost::uniform_on_sphere<double> s(3); // that's my distribution
			boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> > mysphere(rng,s);
	
			Info::Inst()->write(string("Generating orientations for orientational averaging using sphere,boost_uniform_on_sphere"));
    
			for(size_t i = 0; i < resolution; ++i)
			{
				vector<double> r = mysphere();
				this->push_back( CartesianCoor3D(r[0],r[1],r[2]) );							
			}
			// broken:
//		} else if (algorithm=="quad") {
//    		double x0,x1,x2,x3;
//            size_t num = 0;
//            //boost::mt19937 rng; // that's my random number generator			
//            //boost::uniform_real<double> s(0,1); // that's my distribution
//			//boost::variate_generator<boost::mt19937, boost::uniform_real<double> > mygen(rng,s);
//			
//    		while (num<resolution) {
//                x0 = (float)rand()/(float)RAND_MAX;
//                x1 = (float)rand()/(float)RAND_MAX;
//                x2 = (float)rand()/(float)RAND_MAX;
//                x3 = (float)rand()/(float)RAND_MAX;
//                
//    			//x0 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
//    			//x1 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
//    			//x2 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
//    			//x3 = (2.0*(rand()*1.0/RAND_MAX)) - 1.0;
//    			double xl = powf(x0,2) + powf(x1,2) + powf(x2,2) + powf(x3,2);
//    			if ( xl >= 1.0 ) continue;
//
//                double rsign = (float)rand()/(float)RAND_MAX;
//                double sign = 1;
//                if (rsign<0.5) sign=-1;
//                
//    			double x = 2* (x1*x3+x0*x2) / xl;
//    			double y = 2* (x2*x3-x0*x1) / xl;
//    			double z = 2* (powf(x0,2)+powf(x3,2)-powf(x1,2)-powf(x2,2)) / xl;
//
//    			SphericalCoor3D s = CartesianCoor3D(x,y,z);
//    			// take length from scattering vector and angles from unisphere		
//    			this->push_back( CartesianCoor3D(SphericalCoor3D(s.r,s.phi,sign*s.theta)) );				
//    			num++;
//    	    }
		} else {
			Err::Inst()->write(string("Vectors algorithm not understood: ")+algorithm);
			throw;
		}

		
	} else if (type=="cylinder") {
		if (algorithm=="boost_uniform_on_sphere") {
			boost::mt19937 rng; // that's my random number generator
			rng.seed(seed);
			boost::uniform_on_sphere<double> s(2); // that's my distribution
			boost::variate_generator<boost::mt19937, boost::uniform_on_sphere<double> > mysphere(rng,s);
	
			Info::Inst()->write(string("Generating orientations for orientational averaging using cylinder,boost_uniform_on_sphere"));
	
			for(size_t i = 0; i < resolution; ++i)
			{
				vector<double> r = mysphere();
				this->push_back( CartesianCoor3D(r[0],r[1],0) ); 							
			}			
		} else if (algorithm=="raster_linear") {
	
			const double M_2PI = 2*M_PI;
			const double radincr = (M_2PI)/(360*resolution);			

			Info::Inst()->write(string("Generating orientations for orientational averaging using cylinder,raster_linear"));

			for (double phi=0;phi<M_2PI;phi+=radincr) {	
					this->push_back( CartesianCoor3D(cos(phi),sin(phi),0) );	
			}
		} else {
			Err::Inst()->write(string("Vectors algorithm not understood: ")+algorithm);
			throw;
		}
	} else {
		Err::Inst()->write(string("Vectors orientation averaging type not understood: ")+type);
		throw;
	}
	
	Info::Inst()->write(string("Initialized orientational averaging with ")+boost::lexical_cast<string>(size())+string(" vectors."));
}


void ScatteringAverageOrientationMultipoleMomentsParameters::create() {
	
	if (type=="file") {
		Info::Inst()->write(string("Reading multipole moments for orientational averaging from file: ")+filepath);
        
		ifstream mmfile(filepath.c_str());
	
        long major, minor; 
		while (mmfile >> major >> minor) {
			this->push_back(std::make_pair(major,minor));
		}
		
	} else if (type=="resolution") {
		Info::Inst()->write(string("Generating multipole moments for orientational averaging using a maxium major=")+boost::lexical_cast<string>(resolution));

	    if (Params::Inst()->scattering.average.orientation.multipole.type=="sphere") {
            this->push_back(std::make_pair(0L,0L));
            for(long l = 1; l <= resolution; ++l)
            {
                for(long m = -l; m <= l; ++m)
                {
                    this->push_back(std::make_pair(l,m));
                }
            }

	    } else if (Params::Inst()->scattering.average.orientation.multipole.type=="cylinder") {
	        this->push_back(std::make_pair(0L,0L));
            for(long l = 1; l <= resolution; ++l)
            {
                this->push_back(std::make_pair(l,0L));
                this->push_back(std::make_pair(l,1L));
                this->push_back(std::make_pair(l,2L));
                this->push_back(std::make_pair(l,3L));
            }
            
	    } else {
            Err::Inst()->write(string("Type not understood: scattering.average.orientation.multipole.moments.type=")+type);
            throw;
	    }
    } else {
        Err::Inst()->write(string("Type not understood: scattering.average.orientation.multipole.type=")+Params::Inst()->scattering.average.orientation.multipole.type);
        throw;
    }

    // check validity of moments
    Info::Inst()->write(string("Checking multipole moments."));
    if (Params::Inst()->scattering.average.orientation.multipole.type=="cylinder") {
        for(size_t i = 0; i < size(); ++i)
        {
            std::pair<long,long> mm = this->at(i);
            if (mm.first<0) {            
                Err::Inst()->write(string("Major multipole moment must be >= 0!"));
                Err::Inst()->write(string("major=")+boost::lexical_cast<string>(mm.first)+string(", minor=")+boost::lexical_cast<string>(mm.second));                
                throw;                
            }
            if ((mm.second<0) || (mm.second>3) ) {
                Err::Inst()->write(string("Minor multipole moment must be between 0 and 3!"));
                Err::Inst()->write(string("major=")+boost::lexical_cast<string>(mm.first)+string(", minor=")+boost::lexical_cast<string>(mm.second));                
                throw;
            }
            if ((mm.first==0) && (mm.second!=0) ) {
                Err::Inst()->write(string("Minor multipole moment must be 0 for Major 0!"));
                Err::Inst()->write(string("major=")+boost::lexical_cast<string>(mm.first)+string(", minor=")+boost::lexical_cast<string>(mm.second));
                throw;
                throw;
            }
        }
    } else if (Params::Inst()->scattering.average.orientation.multipole.type=="sphere") {
        for(size_t i = 0; i < size(); ++i)
        {
            std::pair<long,long> mm = this->at(i);
            if (mm.first<0) {            
                Err::Inst()->write(string("Major multipole moment must be >= 0!"));
                Err::Inst()->write(string("major=")+boost::lexical_cast<string>(mm.first));
                throw;                
            }            
            if (labs(mm.second)>mm.first) {
                Err::Inst()->write(string("Minor multipole moment must be between -Major and +Major!"));
                Err::Inst()->write(string("major=")+boost::lexical_cast<string>(mm.first)+string(", minor=")+boost::lexical_cast<string>(mm.second));
                throw;
            }
        }
    }
    
	Info::Inst()->write(string("Initialized orientational averaging with ")+boost::lexical_cast<string>(size())+string(" moments."));
}


void ScatteringVectorsParameters::create_from_scans() {
	// read out the scans and push the results onto the internal vector

	if (scans.size()>3) {
		Err::Inst()->write("More than 3 scan definitions are not supported.");
		throw;
	}
	
	// local vector unfolds each scan first, then we'll do element-wise vector addition
	vector< vector<CartesianCoor3D> > qvectors(scans.size());
	for(size_t i = 0; i < scans.size(); ++i)
	{
        if (scans[i].points==0) continue;
        if (scans[i].points==1) {
            double scal = (scans[i].from+scans[i].to)/2;
			qvectors[i].push_back(scal*scans[i].basevector);            
            continue;
        }
        if (scans[i].points==2) {
            qvectors[i].push_back(scans[i].from*scans[i].basevector);            
    		qvectors[i].push_back(scans[i].to*scans[i].basevector);            
            continue;
        }
        // else
        qvectors[i].push_back(scans[i].from*scans[i].basevector);            
        for (size_t j=1;j<(scans[i].points-1);j++) {		
			double scal = scans[i].from + powf((j*1.0/(scans[i].points-1)),scans[i].exponent)*(scans[i].to-scans[i].from);
			qvectors[i].push_back(scal*scans[i].basevector);
		}		
		qvectors[i].push_back(scans[i].to*scans[i].basevector);            
	}

	// trivial case: only one scan!
	if (scans.size()==1) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			this->push_back(qvectors[0][i]);
		}
		return;
	}
	// if scans.size()>1 , not trivial..
	
	if (scans.size()==2) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			for(size_t j = 0; j < qvectors[1].size(); ++j)
			{
				this->push_back(qvectors[0][i]+qvectors[1][j]);
			}
		}	
	}
	
	if (scans.size()==3) {
		for(size_t i = 0; i < qvectors[0].size(); ++i)
		{
			for(size_t j = 0; j < qvectors[1].size(); ++j)
			{
				for(size_t k = 0; k < qvectors[2].size(); ++k)
				{
					this->push_back(qvectors[0][i]+qvectors[1][j]+qvectors[2][k]);
				}
			}
		}	
	}	
}

// end of file
