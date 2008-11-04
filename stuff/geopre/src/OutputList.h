#ifndef OUTPUTLIST_H_
#define OUTPUTLIST_H_

#include "global.h"
#include "OutputFile.h"
#include <vector>
#include "Scene.h"

class Scene;
class OutputFile;

using namespace std;

class OutputList
{
public:
	OutputList();
	virtual ~OutputList();
	
	void generateOutput();

	string sParameter;
	string sGrid;
	string sFormat;
	string sObjectCut;
	bool bYeeGrid_isset;
	bool bYeeGrid;
	bool bUseCluster_isset;
	bool bUseCluster;
	
	vector<OutputFile*> vOutputFiles;
	Scene* ptScene;	 
	string sName;
};

#endif /*OUTPUTLIST_H_*/
