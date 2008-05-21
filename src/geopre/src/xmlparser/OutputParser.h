#ifndef OUTPUTPARSER_H_
#define OUTPUTPARSER_H_

#include <stdio.h>
#include <cstdlib>
#include "expatpp.h"
#include <string.h>
#include <vector>
#include <map>
#include "../Scene.h"
#include "../constants.h" 
#include "../OutputList.h"

using namespace std;

class OutputParser : public expatppNesting
{
public:
	OutputParser(expatppNesting* parent);
	
	virtual ~OutputParser();
		
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);

	Scene* sScene;
	
protected:
	bool m_bDefaultOLSet;
	bool m_bDefaultOFSet;
	string m_sCurrentOL;
};

#endif /*OUTPUTPARSER_H_*/
