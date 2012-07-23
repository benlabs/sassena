/** \file
This file contains an XML Interface class, which re-interpretates the c style libxml2 interface and maps them to C++ elements. It is meant to parse an XML through the use of XPATH.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

// direct header
#include "io/xml_interface.hpp"
#include <assert.h>

using namespace std;




template<> bool XMLInterface::get_value<bool>(const char* xpathexp) {
	std::vector<XMLElement> elements = get(xpathexp);
	// if elements has more than one entry, then xpathexp is ambigious

	if (elements.size()==0) {
		return boost::lexical_cast<bool>("");
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
	std::string result_upper = result;
	boost::to_upper(result_upper);
	// circumvent differences b/w lexical_cast and xml specification
	if (result_upper=="FALSE") result = "0";
	if (result_upper=="TRUE") result = "1";		
	return boost::lexical_cast<bool>(result);
}

XMLElement::XMLElement(xmlNodePtr ptr) {
	node_ptr = ptr;
	if (node_ptr->children!=NULL) { 
		xmlNodePtr cur = node_ptr->children; 
		while(cur!=NULL) { 
			children.push_back(cur);
			cur=cur->next;
		}
	}
}

void XMLElement::print(std::string prepend, bool showchildren) {
	cout << prepend << "type   = " << this->type() << endl;
	cout << prepend << "name   = " << this->name() << endl;
	cout << prepend << "content= " << this->content() << endl;
	if (!showchildren) return;
	if (this->children.size()==0) return;
	
	cout << prepend << "children: [" << this->children.size() << "]" << endl;
	for(size_t i = 0; i < this->children.size(); ++i)
	{
		cout << prepend << "-----" << endl;			
		this->children[i].print(string(prepend)+"  ", showchildren);
		cout << prepend << "-----" << endl;	
	}
}

XMLInterface::XMLInterface(std::string filename) {
	xmlInitParser();
	try { 	
		p_doc = xmlParseFile(filename.c_str()); 
	} 
	catch(...) { 
		Err::Inst()->write(string("Parse error"));
		throw;
	} 
		
	/* Create xpath evaluation context */
    p_xpathCtx = xmlXPathNewContext(p_doc);
    if(p_xpathCtx == NULL) {
		Err::Inst()->write(string("Error: unable to create new XPath context"));
        xmlFreeDoc(p_doc); 
		throw;
    }

    int maxrecursion=20;
    while(xmlXIncludeProcess(p_doc)!=0) {
        maxrecursion--;
        if (maxrecursion<1) {
            Err::Inst()->write("Inclusion depth higher then maxrecursion level(20). Aborting...");
            throw;
        }
    }
//    xmlDocFormatDump(stdout,p_doc,1) ;
    
}

void XMLInterface::dump(std::vector<char>& c) { 
	int chars;
	xmlChar * data;
	xmlDocDumpMemory(p_doc,&data,&chars);
	std::stringstream ss;
	for(int i = 0; i < chars; ++i) c.push_back(data[i]);
	xmlFree(data);
}

XMLInterface::~XMLInterface() { 
	
	xmlXPathFreeContext(p_xpathCtx); 
	xmlFreeDoc(p_doc); 
	xmlCleanupParser(); 
	
}


//inline std::string to_s(const xmlChar* value) {
//	std::stringstream ss;
//	ss << value;
//	return ss.str();	
//}
//
//inline std::string to_s(const xmlNs* value) {
//	std::stringstream ss;
//	ss << value;
//	return ss.str();	
//}


vector<XMLElement> XMLInterface::get(const char* xpathexp) {
	return get(string(xpathexp));	
}

bool XMLInterface::exists(const char* xpathexp) {
	const xmlChar* p_xpath_exp = xmlCharStrdup(xpathexp); 
	xmlXPathObjectPtr p_xpathObj = xmlXPathEvalExpression(p_xpath_exp, p_xpathCtx);
	if (p_xpathObj->nodesetval->nodeNr==0) {
		return false;
	} else {
		return true;
	}
}

vector<XMLElement> XMLInterface::get(std::string xpathexp) { 
	// PREPARE XPATH expressions first...
	const xmlChar* p_xpath_exp = xmlCharStrdup(xpathexp.c_str()); 

	xmlXPathObjectPtr p_xpathObj = xmlXPathEvalExpression(p_xpath_exp, p_xpathCtx);
	if (p_xpathObj==NULL) {
		Err::Inst()->write(string("XPath expression not found: ")+xpathexp);
		throw;
	}
	
	xmlNodeSetPtr nodes = p_xpathObj->nodesetval; // nodes	
	size_t size = (nodes) ? nodes->nodeNr : 0;

	vector<XMLElement> result;	
	for(size_t i = 0; i < size; ++i)
	{
		result.push_back(nodes->nodeTab[i]);
		
	}
    xmlXPathFreeObject(p_xpathObj);

	return result;
}


// end of file
