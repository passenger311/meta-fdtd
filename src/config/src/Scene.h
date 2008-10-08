#ifndef SCENE_H_
#define SCENE_H_

#include <vector>
#include "errorclasses.h"
#include "objects/CObject.h"
#include "OutputFile.h"
#include "OutputList.h"
#include "GridBox.h"
#include <string>
#include <map>

class OutputFile;
class OutputList;
class GridBox;

using namespace std;

class Scene
{
public:
	vector<CObject*> objects;
	map<string, CObject*, stringCmp> moNamedObjects;
	map<string, GridBox*, stringCmp> mgNamedGrids;
	map<string, OutputFile*, stringCmp> moNamedOutputFiles;
	map<string, OutputList*, stringCmp> moNamedOutputLists;

	string sCurrentClusterNode;
	bool bSilentMode;
	bool bBackupOldFiles;
	bool bOverwriteFiles;
	bool bNoErrorOutput;
	bool bFullStatistics;
	bool bQuit;
	bool bSingleOutputFile;
	
	OutputFile* ptCurrentOutputFile;
	OutputList* ptCurrentOutputList;

	GridBox* ptCurrentGrid;
	
	Exception* exFileReading;
	
public:
	Scene();
	virtual ~Scene();

	void generateOutput();
	void printStatistics();
	void silentMode(bool mode); 
	void backupOldFiles(bool mode);
	void overwriteFiles(bool mode);
	void noErrorOutput(bool mode);
	void fullStatistics(bool mode);
	
	void selectGrid(string gridname,bool exceptionOnUnknown = false);
	void selectOutputList(string olname,bool exceptionOnUnknown = false);
	void selectOutputFile(string ofname,bool exceptionOnUnknown = false);
};

#endif /*SCENE_H_*/
