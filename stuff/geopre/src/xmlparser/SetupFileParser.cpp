#include "SetupFileParser.h"
#include "GeometryParser.h"
#include "GridBoxParser.h"
#include "OutputParser.h"
#include "UnitParser.h"
#include "ClusterParser.h"
#include "ParameterParser.h"

SetupFileParser::SetupFileParser(Scene* scene)
{
	sScene = scene;
}
SetupFileParser::~SetupFileParser()
{
}

void SetupFileParser::startElement(const char* name, const char** atts)
{
	if (sScene->exFileReading != NULL)
		return;
	if (strcmp(name,"Geometry") == 0) {
		GeometryParser *par = new GeometryParser(this);
		par->sScene = sScene;
	} else if (strcmp(name,"GridDefinitions") == 0) {
		GridBoxParser *par = new GridBoxParser(this);
		par->sScene = sScene;
	} else if (strcmp(name,"OutputTemplates") == 0) {
		OutputParser *par = new OutputParser(this);
		par->sScene = sScene;
	}
}

void SetupFileParser::endElement(const XML_Char* name)
{
	
}


void SetupFileParser::charData(const XML_Char *s, int len)
{
	const int leadingSpace = skipWhiteSpace(s);
	if (len==0 || len==leadingSpace)
	return;  // called with whitespace between elements
}
