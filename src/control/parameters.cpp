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
#include <boost/random/mersenne_twister.hpp>	
#include <boost/random/uniform_on_sphere.hpp>	
#include <boost/random/variate_generator.hpp>	

// other headers
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

void Params::write_xml(std::string filename) {
    std::ofstream conf(filename.c_str());
    
    conf.close();
}

void Params::read_xml(std::string filename) {
	
	// store the configuration line by line into an internal buffer, 
	// this is for keeping history
	ifstream conffile(filename.c_str());
	string line;
	while (getline(conffile,line)) {
		carboncopy.push_back(line);
	}
	conffile.close();
		
	// START OF sample section //	
	XMLInterface xmli(filename);
    
	Info::Inst()->write("Reading parameters...");
	
	// now read the parameters
	
	if (xmli.exists("//sample")) {
    	if (xmli.exists("//sample/structure")) {
    	    if (xmli.exists("//sample/structure/file")) {
            	sample.structure.file   = get_filepath(xmli.get_value<std::string>("//sample/structure/file"));
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
                        xmli.set_current(indexes[i]);
                        ids.push_back(xmli.get_value<size_t>("."));
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
                    string filename = "selection.pdb";
                    string format = "pdb";
                    string selector = "beta";
                    string expression = "1|1\\.0|1\\.00";
                    
                    if (xmli.exists("./format")) format = xmli.get_value<string>("./format") ;

                    // this is a convenience overwrite for the ndx file format
                    if (format=="ndx") {
                        filename = "index.ndx";
                        selector = "name";
                        expression = ".*";
                    }
                    if (xmli.exists("./file")) filename = xmli.get_value<string>("./file") ;
                    
                    if (xmli.exists("./selector")) selector = xmli.get_value<string>("./selector") ;
                    if (xmli.exists("./expression")) expression = xmli.get_value<string>("./expression") ;
                    
                    sample.selections[sname] = new SampleFileSelectionParameters(filename,format,selector,expression);                    
                    Info::Inst()->write(string("Creating File Atomselection: ")+sname+string(" , file: ")+filename+string(", format: ")+format+string(", selector: ")+selector+string(", expression: ")+expression);
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
                fset.type = def_format;
                fset.clones = def_clones;
        		fset.filename = get_filepath(xmli.get_value<string>("./file"));
                boost::filesystem::path index_path = get_filepath(xmli.get_value<string>("./file"));
                fset.index = index_path.parent_path().string() +string("/")+ index_path.stem() + (".tnx");
        		if (xmli.exists("./format"))  fset.type = xmli.get_value<string>("./format");
        		if (xmli.exists("./first"))   fset.first  = xmli.get_value<size_t>("./first");
        		if (xmli.exists("./last"))  { fset.last   = xmli.get_value<size_t>("./last"); fset.last_set = true; }
        		if (xmli.exists("./stride"))  fset.stride = xmli.get_value<size_t>("./stride");
        		if (xmli.exists("./clones"))  fset.clones = xmli.get_value<size_t>("./clones");
        		if (xmli.exists("./index"))  fset.index = get_filepath(xmli.get_value<std::string>("./index"));

        		sample.framesets.push_back(fset);
        		Info::Inst()->write(string("Added frames from ")+fset.filename+string(" using format: ")+fset.type);
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
	    		if (xmli.exists("./type"))   motion.type  = xmli.get_value<string>("./type");
	    		if (xmli.exists("./displace"))  motion.displace   = xmli.get_value<double>("./displace");
	    		if (xmli.exists("./frequency"))  motion.frequency   = xmli.get_value<double>("./frequency");			
	    		if (xmli.exists("./seed"))  motion.seed   = xmli.get_value<long>("./seed");			
	    		if (xmli.exists("./sampling"))  motion.seed   = xmli.get_value<long>("./sampling");			
	    		if (xmli.exists("./selection"))  motion.selection   = xmli.get_value<string>("./selection");			
	    		if (xmli.exists("./direction")) {
	    			if (xmli.exists("./direction/x")) motion.direction.x   = xmli.get_value<double>("./direction/x");
	    			if (xmli.exists("./direction/y")) motion.direction.y   = xmli.get_value<double>("./direction/y");
	    			if (xmli.exists("./direction/z")) motion.direction.z   = xmli.get_value<double>("./direction/z");				
	    		} 
                motion.radius = motion.displace*10;
	    		if (xmli.exists("./radius"))  motion.radius   = xmli.get_value<double>("./radius");			
        
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

    	    	sample.alignments.push_back(alignment);
    	    	Info::Inst()->write(string("Adding additional alignment to sample: type=")+alignment.type+string(", selection=")+alignment.selection+string(", order=")+alignment.order);
    	    }
        }	
	
	}	
	
	scattering.type="all";
    scattering.target = "system";
	scattering.background.factor = 0.0;

    scattering.dsp.type="autocorrelate";
    scattering.dsp.method="fftw";
	// defaults
	scattering.average.orientation.type = "none";
	scattering.average.orientation.vectors.type = "sphere";
	scattering.average.orientation.vectors.algorithm = "boost_uniform_on_sphere";
	scattering.average.orientation.vectors.resolution = 100;
	scattering.average.orientation.vectors.seed = 0;
	scattering.average.orientation.vectors.axis = CartesianCoor3D(0,0,1);
	scattering.average.orientation.vectors.file = "qvector-orientations.txt";
	scattering.average.orientation.multipole.type = "sphere";
	scattering.average.orientation.multipole.resolution = 20;
	scattering.average.orientation.multipole.axis = CartesianCoor3D(0,0,1);
	scattering.average.orientation.exact.type = "sphere";

    scattering.signal.file = "signal.h5";
    scattering.signal.fqt = true;
    scattering.signal.fq0 = true;
    scattering.signal.fq = true;
    scattering.signal.fq2 = true;

    if (xmli.exists("//scattering")) {
    	if (xmli.exists("//scattering/type")) {
    		scattering.type = xmli.get_value<string>("//scattering/type");
    	}
    	Info::Inst()->write(string("scattering.type=")+scattering.type);


    	if (xmli.exists("//scattering/target")) {
    		scattering.target = xmli.get_value<string>("//scattering/target");
    	}
    	Info::Inst()->write(string("scattering.target=")+scattering.target);


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
	    				scattering.average.orientation.vectors.seed = xmli.get_value<long>("//scattering/average/orientation/vectors/seed");
                        Info::Inst()->write(string("scattering.average.orientation.vectors.seed=")+boost::lexical_cast<string>(scattering.average.orientation.vectors.seed));
	    			}	
	    			if (xmli.exists("//scattering/average/orientation/vectors/file")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.file = get_filepath(xmli.get_value<string>("//scattering/average/orientation/vectors/file"));
	    			}		
	    			if (xmli.exists("//scattering/average/orientation/vectors/axis")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.vectors.axis.x = xmli.get_value<double>("//scattering/average/orientation/vectors/axis/x");
	    				scattering.average.orientation.vectors.axis.y = xmli.get_value<double>("//scattering/average/orientation/vectors/axis/y");
	    				scattering.average.orientation.vectors.axis.z = xmli.get_value<double>("//scattering/average/orientation/vectors/axis/z");					
	    			}						
	    			
	    			scattering.average.orientation.vectors.create();
	    			
	    		} else if (scattering.average.orientation.type=="multipole") {			    
	    			if (xmli.exists("//scattering/average/orientation/multipole/type")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.multipole.type = xmli.get_value<string>("//scattering/average/orientation/multipole/type");
                        Info::Inst()->write(string("scattering.average.orientation.multipole.type=")+scattering.average.orientation.multipole.type);				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/multipole/resolution")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.multipole.resolution = xmli.get_value<long>("//scattering/average/orientation/multipole/resolution");
                        Info::Inst()->write(string("scattering.average.orientation.multipole.resolution=")+boost::lexical_cast<string>(scattering.average.orientation.multipole.resolution));				    
	    			}
	    			if (xmli.exists("//scattering/average/orientation/multipole/axis")) { // count vectors ... , or order for multipole...
	    				scattering.average.orientation.multipole.axis.x = xmli.get_value<double>("//scattering/average/orientation/multipole/axis/x");
	    				scattering.average.orientation.multipole.axis.y = xmli.get_value<double>("//scattering/average/orientation/multipole/axis/y");
	    				scattering.average.orientation.multipole.axis.z = xmli.get_value<double>("//scattering/average/orientation/multipole/axis/z");					
	    			}						
	    		} else if (scattering.average.orientation.type!="none") {
	    			Err::Inst()->write(string("Orientation averaging type not understood: ")+scattering.average.orientation.type);
	    			throw;
	    		}
	    	}
	    }
	    
	    if (xmli.exists("//scattering/signal")) {
	        if (xmli.exists("//scattering/signal/file")) {
        		scattering.signal.file = xmli.get_value<string>("//scattering/signal/file");	        
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
	    Info::Inst()->write(string("scattering.signal.fqt=")+boost::lexical_cast<string>(scattering.signal.fqt));	
	    Info::Inst()->write(string("scattering.signal.fq0=")+boost::lexical_cast<string>(scattering.signal.fq0));	            
	    Info::Inst()->write(string("scattering.signal.fq=")+boost::lexical_cast<string>(scattering.signal.fq));	
	    Info::Inst()->write(string("scattering.signal.fq2=")+boost::lexical_cast<string>(scattering.signal.fq2));	    
    }  
	// END OF scattering section //

	// START OF limits section //

    // assign default memory limits:
    limits.stage.nodes = 1;
    limits.stage.mode = "mod";
    limits.stage.sync_barrier = 1000;
    limits.stage.memory.buffer = 100*1024*1024;
    limits.stage.memory.data   = 500*1024*1024;
    
    limits.signal.chunksize = 10000;
        
    limits.computation.threads.on = false;
    limits.computation.threads.scatter = 1;
    limits.computation.threads.scatter_timeout = 25; // 25 ms
    limits.computation.threads.dsp = 1;
    limits.computation.memory.signal = 200*1024*1024; // 100MB 

    limits.services.signal.memory.server = 100*1024*1024; // 100MB
    limits.services.signal.memory.client = 10*1024*1024; // 10MB
    limits.services.signal.times.serverflush = 600; // 600 seconds
    limits.services.signal.times.clientflush = 600; // 600 seconds
    
    limits.decomposition.utilization = 0.95; // 5% max loss
    limits.decomposition.partitions.automatic = true; // pick number of independent partitions based on some heuristics
    limits.decomposition.partitions.size = 1; // not used if automatic = true, if false -> this determines the partition size

	if (xmli.exists("//limits")) {       	
	    if (xmli.exists("//limits/stage")) {
        	if (xmli.exists("//limits/stage/nodes")) {
    	        limits.stage.nodes = xmli.get_value<size_t>("//limits/stage/nodes");
	        }	        
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
        	if (xmli.exists("//limits/computation/threads")) {
            	if (xmli.exists("//limits/computation/threads/on")) {
        	        limits.computation.threads.on = xmli.get_value<bool>("//limits/computation/threads/on");
    	        }
            	if (xmli.exists("//limits/computation/threads/scatter")) {
        	        limits.computation.threads.scatter = xmli.get_value<size_t>("//limits/computation/threads/scatter");
    	        }
            	if (xmli.exists("//limits/computation/threads/dsp")) {
        	        limits.computation.threads.dsp = xmli.get_value<size_t>("//limits/computation/threads/dsp");
    	        }
            	if (xmli.exists("//limits/computation/threads/scatter_timeout")) {
        	        limits.computation.threads.scatter_timeout = xmli.get_value<size_t>("//limits/computation/threads/scatter_timeout");
    	        }
	        }
        	if (xmli.exists("//limits/computation/memory")) {
            	if (xmli.exists("//limits/computation/memory/signal")) {
    	            limits.computation.memory.signal = xmli.get_value<size_t>("//limits/computation/memory/signal");
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
    debug.iowrite.buffer = false;
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
};

void Params::init(std::string filename) {
	
	// make the directory of the main configuration file the root for all others
	if (boost::filesystem::path(filename).is_complete()) 
		config_rootpath = boost::filesystem::path(filename).parent_path().string();
	else 
		config_rootpath = ( boost::filesystem::initial_path() / boost::filesystem::path(filename).parent_path() ).string();		

	Info::Inst()->write(string("Looking for configuration file: ") + filename);
	
	read_xml(filename);
}

SampleParameters::~SampleParameters() {
    for (std::map<std::string,SampleSelectionParameters*>::iterator i=selections.begin();i!=selections.end();i++) {
        delete i->second;
    }
}

void ScatteringAverageOrientationVectorsParameters::create() {
	
	if (type=="file") {
		
		ifstream qqqfile(file.c_str());
	
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
	
			for(size_t i = 0; i < resolution; ++i)
			{
				vector<double> r = mysphere();
				this->push_back( CartesianCoor3D(r[0],r[1],r[2]) );							
			}
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
	
			for(size_t i = 0; i < resolution; ++i)
			{
				vector<double> r = mysphere();
				this->push_back( CartesianCoor3D(r[0],r[1],0) ); 							
			}			
		} else if (algorithm=="raster_linear") {
	
			const double M_2PI = 2*M_PI;
			const double radincr = (M_2PI)/(360*resolution);			

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
