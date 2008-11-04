#ifndef SETUPFILEPARSER_H_
#define SETUPFILEPARSER_H_

#include "expatpp.h"

#include "../Scene.h"

class Scene;

class SetupFileParser : public expatppNesting
{
public:
	Scene* sScene;
public:
	SetupFileParser(Scene* scene);
	virtual ~SetupFileParser();

	virtual void startElement(const XML_Char *name, const XML_Char **atts);
	virtual void endElement(const XML_Char* name);
	virtual void charData(const XML_Char *s, int len);
};

#endif /*SETUPFILEPARSER_H_*/
