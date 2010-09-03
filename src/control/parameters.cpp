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
	
	    
	sample.structure.file   = get_filepath(xmli.get_value<std::string>("//sample/structure/file"));
	sample.structure.format = xmli.get_value<std::string>("//sample/structure/format");
	vector<XMLElement> selections = xmli.get("//sample/selections/selection");
	for(size_t i = 0; i < selections.size(); ++i)
	{
		xmli.set_current(selections[i]);
		double sv = 1.0;
        string sn = "beta";
        string ff = "pdb";
        string fn = "selection";
        string sname = string("selection_")+boost::lexical_cast<string>(i);
		if (xmli.exists("./name")) {
			sname = xmli.get_value<string>("./name") ;
		}
		if (xmli.exists("./file")) {		
			fn = get_filepath( xmli.get_value<string>("./file") );
		}
		if (xmli.exists("./format")) {		
			ff = xmli.get_value<string>("./format") ;
		}
		if (xmli.exists("./select")) {		
			sn = xmli.get_value<string>("./select") ;
		}
		if (xmli.exists("./select_value")) {
			sv = xmli.get_value<double>("./select_value") ;
		}
		sample.groups[sname] = SampleGroupParameters(sname,fn,ff,sn,sv);
	}
	
	// END OF sample section //	

	// START OF scattering section //	

	// END OF scattering section //	

	// START OF output section //	

	// END OF output section //	

	// START OF limits section //	

	// END OF limits section //	

	// START OF debug section //	

	// END OF debug section //	
	
	// START OF sample section //	

	// read in frame information
	size_t def_first=0;	
	size_t def_last=0;	
	bool def_last_set=false; 
	size_t def_stride = 1;
    size_t def_clones = 1;
	if (xmli.exists("//sample/frames/first"))   def_first  = xmli.get_value<size_t>("//sample/frames/first");
	if (xmli.exists("//sample/frames/last"))  { def_last   = xmli.get_value<size_t>("//sample/frames/last"); def_last_set = true; }
	if (xmli.exists("//sample/frames/stride"))  def_stride = xmli.get_value<size_t>("//sample/frames/stride");
	if (xmli.exists("//sample/frames/clones"))  def_clones = xmli.get_value<size_t>("//sample/frames/clones");
	
	vector<XMLElement> framesets = xmli.get("//sample/frames/frameset");
	for(size_t i = 0; i < framesets.size(); ++i)
	{
		xmli.set_current(framesets[i]);
		SampleFramesetParameters fset;	
		fset.first = def_first;	fset.last = def_last; fset.last_set = def_last_set; fset.stride = def_stride;
        fset.clones = def_clones;
		fset.filename = get_filepath(xmli.get_value<string>("./file"));
		fset.type = xmli.get_value<string>("./format");
		if (xmli.exists("./first"))   fset.first  = xmli.get_value<size_t>("./first");
		if (xmli.exists("./last"))  { fset.last   = xmli.get_value<size_t>("./last"); fset.last_set = true; }
		if (xmli.exists("./stride"))  fset.stride = xmli.get_value<size_t>("./stride");
		if (xmli.exists("./clones"))  fset.clones = xmli.get_value<size_t>("./clones");
		
		sample.frames.push_back(fset);
		Info::Inst()->write(string("Added frames from ")+fset.filename+string(" using format: ")+fset.type);
		Info::Inst()->write(string("Options: first=")+to_s(fset.first)+string(", last=")+to_s(fset.last)+string(", lastset=")+to_s(fset.last_set)+string(", stride=")+to_s(fset.stride));		
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
			motion.frequency=2*M_PI/1000.0; // corresponds to one full cycle per 1000 frames, used for linear oscillation and rotation
			if (xmli.exists("./type"))   motion.type  = xmli.get_value<string>("./type");
			if (xmli.exists("./displace"))  motion.displace   = xmli.get_value<double>("./displace");
			if (xmli.exists("./frequency"))  motion.frequency   = xmli.get_value<double>("./frequency");			
			if (xmli.exists("./seed"))  motion.seed   = xmli.get_value<long>("./seed");			
			if (xmli.exists("./sampling"))  motion.seed   = xmli.get_value<long>("./sampling");			
			if (xmli.exists("./selection"))  motion.selection   = xmli.get_value<string>("./selection");			
			if (xmli.exists("./direction")) {
				motion.direction.x   = xmli.get_value<double>("./direction/x");
				motion.direction.y   = xmli.get_value<double>("./direction/y");
				motion.direction.z   = xmli.get_value<double>("./direction/z");				
			} 

			sample.motions.push_back(motion);
			Info::Inst()->write(string("Adding additional motion to sample: type=")+motion.type+string(", displacement=")+to_s(motion.displace)+string(", selection=")+motion.selection);
		}
	}	
		
    
    if (xmli.exists("//sample/alignments")) {

	    vector<XMLElement> alignments = xmli.get("//sample/alignments/alignment");
	    for(size_t i = 0; i < alignments.size(); ++i)
	    {
	    	xmli.set_current(alignments[i]);
	    	SampleAlignmentParameters alignment;	
	    	alignment.type = "center";	 
	    	alignment.selection = "";
            alignment.order = "pre";
	    	if (xmli.exists("./type"))   alignment.type  = xmli.get_value<string>("./type");
	    	if (xmli.exists("./selection"))  alignment.selection   = xmli.get_value<string>("./selection");			
	    	if (xmli.exists("./order"))  alignment.order   = xmli.get_value<string>("./order");			
        
	    	sample.alignments.push_back(alignment);
            string selection_string = "system";
            if (alignment.selection!="") selection_string = alignment.selection;
	    	Info::Inst()->write(string("Adding additional alignment to sample: type=")+alignment.type+string(", selection=")+selection_string+string(", order=")+alignment.order);
	    }
    }	

	// END OF sample section //
	// START OF scattering section //

    scattering.center = false;
	if (xmli.exists("//scattering/center")) {
        scattering.center = xmli.get_value<bool>("//scattering/center");	    
	    Info::Inst()->write(string("scattering.center=")+to_s(scattering.center));		    
	} 
	scattering.background.factor = 0.0;
    	
	if (xmli.exists("//scattering/background")) {
		if (xmli.exists("//scattering/background/factor")) {
		    scattering.background.factor = xmli.get_value<double>("//scattering/background/factor");
		    Info::Inst()->write(string("scattering.background.factor=")+to_s(scattering.background.factor));		    		    
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
    	    	Info::Inst()->write(string("Rescaling volumes: selection=")+kappa.selection+string(", value=")+to_s(kappa.value));
    	    }
        }
	}	
	// generating qqqvectors, i.e. the spectrum

	string vt = xmli.get_value<string>("//scattering/vectors/type");
	if (vt=="single") {	

		double x = xmli.get_value<double>("//scattering/vectors/single/x");
		double y = xmli.get_value<double>("//scattering/vectors/single/y");
		double z = xmli.get_value<double>("//scattering/vectors/single/z");

		scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
	}
	else if (vt=="scan") {	
		
		vector<XMLElement> scans = xmli.get("//scattering/vectors/scan");

		for(size_t i = 0; i < scans.size(); ++i)
		{
			xmli.set_current(scans[i]);
			ScatteringVectorsScanParameters sc;
			
			sc.basevector.x = xmli.get_value<double>("./basevector/x");
			sc.basevector.y = xmli.get_value<double>("./basevector/y");
			sc.basevector.z = xmli.get_value<double>("./basevector/z");

			sc.from     = xmli.get_value<double>("./from");
			sc.to       = xmli.get_value<double>("./to");
			sc.points   = xmli.get_value<size_t>("./points");
			if (xmli.exists("./exponent")) sc.exponent = xmli.get_value<double>("./exponent");
			scattering.qvectors.scans.push_back(sc);
		}
		
		scattering.qvectors.create_from_scans();
	}
	else if (vt=="file") {
		string qqqfilename = get_filepath(xmli.get_value<string>("//scattering/vectors/file"));
		ifstream qqqfile(qqqfilename.c_str());
		
		double x,y,z; 
		while (qqqfile >> x >> y >> z) {
			scattering.qvectors.push_back(CartesianCoor3D(x,y,z));
		}
	}
	
    scattering.correlation.type="instant-time";
    scattering.correlation.method="direct";
    
	if (xmli.exists("//scattering/correlation")) {
		if (xmli.exists("//scattering/correlation/type")) {
			scattering.correlation.type = xmli.get_value<string>("//scattering/correlation/type");
		    Info::Inst()->write(string("scattering.correlation.type=")+scattering.correlation.type);
		}
		if (xmli.exists("//scattering/correlation/method")) {
			scattering.correlation.method = xmli.get_value<string>("//scattering/correlation/method");
		    Info::Inst()->write(string("scattering.correlation.method=")+scattering.correlation.method);
		}
	}

	// defaults
	scattering.average.orientation.type = "none";
	scattering.average.orientation.vectors.type = "sphere";
	scattering.average.orientation.vectors.algorithm = "boost_uniform_on_sphere";
	scattering.average.orientation.vectors.resolution = 100;
	scattering.average.orientation.vectors.seed = 0;
	scattering.average.orientation.vectors.axis = CartesianCoor3D(0,0,1);
	scattering.average.orientation.vectors.file = "";

	scattering.average.orientation.multipole.type = "sphere";
	scattering.average.orientation.multipole.resolution = 20;
	scattering.average.orientation.multipole.axis = CartesianCoor3D(0,0,1);

	scattering.average.orientation.exact.type = "sphere";
		
	if (xmli.exists("//scattering/average")) {
		if (xmli.exists("//scattering/average/orientation")) {
			if (xmli.exists("//scattering/average/orientation/type")) { // vectors multipole
				scattering.average.orientation.type = xmli.get_value<string>("//scattering/average/orientation/type");
			}
		    Info::Inst()->write(string("scattering.average.orientation.type=")+scattering.average.orientation.type);
			
			if (
			    (scattering.average.orientation.type=="vectors") ||
			    (scattering.average.orientation.type=="vectorsthread")
			    ) {
        
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
                    Info::Inst()->write(string("scattering.average.orientation.vectors.resolution=")+to_s(scattering.average.orientation.vectors.resolution));				    
				}
				if (xmli.exists("//scattering/average/orientation/vectors/seed")) { // count vectors ... , or order for multipole...
					scattering.average.orientation.vectors.seed = xmli.get_value<long>("//scattering/average/orientation/vectors/seed");
                    Info::Inst()->write(string("scattering.average.orientation.vectors.seed=")+to_s(scattering.average.orientation.vectors.seed));
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
                    Info::Inst()->write(string("scattering.average.orientation.multipole.resolution=")+to_s(scattering.average.orientation.multipole.resolution));				    
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
	
    scattering.interference.type="all";
	
	if (xmli.exists("//scattering/interference")) {
		if (xmli.exists("//scattering/interference/type")) {
			scattering.interference.type = xmli.get_value<string>("//scattering/interference/type");
		}
	}
	Info::Inst()->write(string("scattering.interference.type=")+scattering.interference.type);

	
    scattering.target = "system";
		
	if (xmli.exists("//scattering/target")) {
		scattering.target = xmli.get_value<string>("//scattering/target");
	}
	Info::Inst()->write(string("scattering.target=")+scattering.target);
	
	// END OF scattering section //
	// START OF output section //
	
	if (!xmli.exists("//output")) {
		Err::Inst()->write("You have to specify an output section");
		throw;
	}
	
	if (xmli.exists("//output/prefix")) {
		output.prefix = xmli.get_value<string>("//output/prefix");
	}
	
	vector<XMLElement> outputfiles = xmli.get("//output/file");
	
	for(size_t i = 0; i < outputfiles.size(); ++i)
	{
		xmli.set_current(outputfiles[i]);
		OutputFileParameters ofp;
		ofp.name = xmli.get_value<string>("./name");
		ofp.method = xmli.get_value<string>("./method");
		ofp.format = xmli.get_value<string>("./format");
		ofp.filename = get_filepath(output.prefix + ofp.name + string(".") + ofp.format  );
		output.files.push_back(ofp);
	}
	
	if (output.files.size()==0) {
		Err::Inst()->write("Aborting. No output files defined.");
		throw;
	}

	// END OF output section //
	// START OF limits section //

    // assign default memory limits:
    limits.memory.scattering_matrix = 100*1024*1024; // 100MB
    limits.memory.coordinate_sets = 500*1024*1024; // 500MB
    limits.decomposition.static_imbalance = 0.05; // 5% max
    limits.decomposition.partitions.automatic = true; // pick number of independent partitions based on some heuristics
    limits.decomposition.partitions.max = 1000000; // virtually dont limit the maximum number of partitions
    limits.decomposition.partitions.count = 1; // not used if automatic = true, if false -> this determines the split factor

	if (xmli.exists("//limits")) {        
    	if (xmli.exists("//limits/memory")) {
        	if (xmli.exists("//limits/memory/scattering_matrix")) {
    	        limits.memory.scattering_matrix = xmli.get_value<size_t>("//limits/memory/scattering_matrix");
	        }
        	if (xmli.exists("//limits/memory/coordinate_sets")) {	        
			    limits.memory.coordinate_sets = xmli.get_value<size_t>("//limits/memory/coordinate_sets");
    	    }
	    }
    	if (xmli.exists("//limits/decomposition")) {
        	if (xmli.exists("//limits/decomposition/static_imbalance")) {
			    limits.decomposition.static_imbalance = xmli.get_value<double>("//limits/decomposition/static_imbalance");
            }    	    
        	if (xmli.exists("//limits/decomposition/partitions")) {
            	if (xmli.exists("//limits/decomposition/partitions/automatic")) {
    			    limits.decomposition.partitions.automatic = xmli.get_value<bool>("//limits/decomposition/partitions/automatic");
                }
            	if (xmli.exists("//limits/decomposition/partitions/max")) {
    			    limits.decomposition.partitions.max = xmli.get_value<size_t>("//limits/decomposition/partitions/max");
                }
            	if (xmli.exists("//limits/decomposition/partitions/count")) {
    			    limits.decomposition.partitions.count = xmli.get_value<size_t>("//limits/decomposition/partitions/count");                    
                }
            }
        }		
	}

	// END OF limits section //
	// START OF debug section //

	
	debug.timer = false; // this adds a log message when a timer is started/stopped
	debug.barriers = false; // this de-/activates collective barriers before each collective operation, this way all nodes are synchronized before the communication takes place. This is an important step towards analysis of timing.
	if (xmli.exists("//debug")) {
		if (xmli.exists("//debug/timer")) {
			debug.timer = xmli.get_value<bool>("//debug/timer");
		}
		if (xmli.exists("//debug/barriers")) {
			debug.barriers = xmli.get_value<bool>("//debug/barriers");
		}
	}
	
	// initialize some of the runtime parameters:
    runtime.limits.cache.coordinate_sets = 1;

};

string Params::guessformat(string filename) {
	// do the best you can to guess the format
	// guess by filename
	boost::regex e_xml(".*\\.xml$");
	
	if (boost::regex_match(filename,e_xml)) return "xml";
	
	// else: problem
	Err::Inst()->write("Parameter file format could not be detected (ending with: xml ?)");
	throw;
}

bool Params::check() {
	return true;
	// implement a sanity check for the parameters
}

void Params::init(std::string filename) {
	string format = guessformat(filename);
	
	// make the directory of the main configuration file the root for all others
	if (boost::filesystem::path(filename).is_complete()) 
		config_rootpath = boost::filesystem::path(filename).parent_path().string();
	else 
		config_rootpath = ( boost::filesystem::initial_path() / boost::filesystem::path(filename).parent_path() ).string();		

	Info::Inst()->write(string("Looking for configuration file: ") + filename);
	
	read_xml(filename);

	Info::Inst()->write("Checking parameters...");	
	
	if (check()) {
		Info::Inst()->write("Check succeeded. Parameters seem OK.");			
	}
	else {
		Err::Inst()->write("Check failed. Please check your input file.");
		throw;		
	}
	
}

void Params::write(std::string filename, std::string format) {

	if (format=="xml") {
		Err::Inst()->write("XML parameters format (write) not yet implemented.");
		throw;
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
