#ifndef GEOMETRYPARSER_H_
#define GEOMETRYPARSER_H_

#include <stdio.h>
#include <cstdlib>
#include "expatpp.h"
#include <string.h>
#include <vector>
#include <map>
#include "../Scene.h"
#include "../geoclasses/algebra3.h" 
#include "../geoclasses/allobjects.h"
#include "../constants.h" 

using namespace std;


class GeometryParser : public expatppNesting {
public:
	GeometryParser()
	{
	}
	GeometryParser(expatppNesting* parent) : expatppNesting(parent) 
	{
	}
	
	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);

	Scene* sScene;

protected:
	vector<CObject*> movCurrentObjects;
	map<const char*, CObject*, strCmp> momObjects;
//	CObject* getCObject(const char* name, dictmap attributes);

};

#endif /*GEOMETRYPARSER_H_*/
