/** \file
This file contains an XML Interface class, which re-interpretates the c style libxml2 interface and maps them to C++ elements. It is meant to parse an XML through the use of XPATH.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef IO__XML_INTERFACE_HPP_
#define IO__XML_INTERFACE_HPP_

// common header
#include "common.hpp"

// standard header
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

// special library headers
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xinclude.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
// other headers
#include "log/log.hpp"

/** 
 Models a XML node and provides convenience functions to access its properties.
*/
class XMLElement {
	xmlNodePtr node_ptr;
public:
	std::vector<XMLElement> children;

	XMLElement(xmlNodePtr ptr);

	void set(xmlNodePtr ptr) { node_ptr = ptr; }

	xmlElementType type()    { return node_ptr->type; }
	std::string    name()    { if (node_ptr->name==NULL) return ""; else return (const char*) node_ptr->name; }
	std::string    content() { if (node_ptr->content==NULL) return ""; else return (const char*)  node_ptr->content; }
	
	void print(std::string prepend, bool showchildren);
	
	xmlNodePtr get_node_ptr() { return node_ptr; }
};

/** 
 Models a XML file and allows access through XPATH.
*/
class XMLInterface {
	
	xmlDocPtr p_doc;
	xmlXPathContextPtr p_xpathCtx; 
	
public:
	XMLInterface(std::string filename);
	~XMLInterface();

	std::vector<XMLElement> get(const char* xpathexp);	
	std::vector<XMLElement> get(std::string xpathexp);

	void dump(std::vector<char>& c);
	
	template<class convT> convT get_value(const char* xpathexp);

	bool exists(const char* xpathexp);
	void set_current(XMLElement thisel) {	p_xpathCtx->node = thisel.get_node_ptr();}
	
	
};


/** 
 Provides type safe access to content element in XML files through XPATH.
*/
template<class convT> convT XMLInterface::get_value(const char* xpathexp) {
	std::vector<XMLElement> elements = get(xpathexp);
	// if elements has more than one entry, then xpathexp is ambigious

	if (elements.size()==0) {
		return boost::lexical_cast<convT>("");
	}
	else if (elements.size()>1) {
		Err::Inst()->write(std::string("Xpathexp is ambigious, multiple fields matched: ")+xpathexp);
		throw;
	}

	XMLElement& thisel = elements[0];
	
	if (thisel.type()!=XML_ELEMENT_NODE) {
		Err::Inst()->write(std::string("Xpathexp doesn't resolve to XML_ELEMENT_NODE: ")+xpathexp);
		throw;
	}

	size_t textelements = 0;
	std::string result;
	for(size_t i = 0; i < thisel.children.size(); ++i)
	{
		if (thisel.children[0].type()==XML_TEXT_NODE) {
			textelements++;
			result = thisel.children[0].content();
		}
	}
	
	// textelements == 0 corresponds to an empty text field. That is legal!
	
	if (textelements>1) {
		Err::Inst()->write(std::string("Xpathexp resolves to more than one text field: ")+xpathexp);
		throw;
	}
	
	// trim whitespaces!
	boost::trim(result);
	// circumvent differences b/w lexical_cast and xml specification

	return boost::lexical_cast<convT>(result);
}


#endif

// end of file
