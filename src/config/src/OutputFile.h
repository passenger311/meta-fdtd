#ifndef OUTPUTFILE_H_
#define OUTPUTFILE_H_

#include "global.h"
#include "Scene.h"

class Scene;

class OutputFile
{
public:
	OutputFile();
	virtual ~OutputFile();
	
	void generateOutput();

public:
	string sParameter;
	string sGrid;
	string sFormat;
	string sFile;
	string sObjectCut;
	bool bYeeGrid;
	bool bUseCluster;
	
	Scene* ptScene;	 
	string sName;
};

#endif /*OUTPUTFILE_H_*/
