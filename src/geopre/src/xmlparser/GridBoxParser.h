#ifndef GRIDBOXPARSER_H_
#define GRIDBOXPARSER_H_

#include "expatpp.h"
#include <stdio.h>
#include <cstdlib>
#include <map>
#include "../Scene.h"
#include "../geoclasses/algebra3.h" 
#include "../constants.h" 

class GridBoxParser : public expatppNesting
{
public:
	GridBoxParser(expatppNesting* parent);
	virtual ~GridBoxParser();
	
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);

	Scene* sScene;

protected:
	vector<CObject*> movCurrentObjects;
	map<const char*, CObject*, strCmp> momObjects;
	bool m_bDefaultGridSet;
};

#endif /*GRIDBOXPARSER_H_*/
