/*
 *  settings.hpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

#ifndef SETTINGS_HPP_
#define SETTINGS_HPP_

// common header
#include "common.hpp"

// standard header
#include <map>
#include <string>
#include <vector>

// special library headers
#include <libconfig.h++>

// other headers
#include "coor3d.hpp"

class Settings {
	static std::string config_rootpath;
	static std::vector<std::string> conf_files;
	static bool conf_files_set;
	
	static std::map<std::string,libconfig::Setting*> settings;
	
	static libconfig::Setting* read_from_configuration(std::string filename,bool expandfilename);
	static void unfold_settings(libconfig::Setting* setting);
public:
	static bool read(int& argc,char** argv);
	
	static libconfig::Setting& get(std::string setting) { return *(settings[setting]); }
	static std::string get_filepath(std::string fname);	
	static double get_volume(std::string atomname);
	static void   get_qqqvectors(std::vector<CartesianCoor3D>& qqqvectors);
};

// this is the stem for wrapping the settings
class RuntimeSettings {
	
};


#endif
