#include "Scene.h"
#include "GridBox.h"
#include "constants.h"
#include "errorclasses.h"
#include <sstream>

Scene::Scene()
{
	sCurrentClusterNode = "";
	bSilentMode = false;
	bBackupOldFiles = false;
	bOverwriteFiles = false;
	bNoErrorOutput = false;
	bFullStatistics = false;
	bQuit = false;
	ptCurrentOutputFile = NULL;
	ptCurrentOutputList = NULL;
	ptCurrentGrid = NULL;
	selectGrid("default");
	selectOutputList("default");
	bSingleOutputFile = false;
	exFileReading = NULL;
}

Scene::~Scene()
{
}

void Scene::silentMode(bool mode) 
{
  bSilentMode = mode;
}

void Scene::backupOldFiles(bool mode)
{
  bBackupOldFiles = mode;
}

void Scene::overwriteFiles(bool mode)
{
  bOverwriteFiles = mode;
}

void Scene::noErrorOutput(bool mode)
{
  bNoErrorOutput = mode;
}

void Scene::fullStatistics(bool mode)
{
  bFullStatistics = mode;
}


void Scene::generateOutput()
{
	if (bSingleOutputFile)
		ptCurrentOutputFile->generateOutput();
	else
		ptCurrentOutputList->generateOutput();
}

void Scene::printStatistics()
{
	
}

void Scene::selectGrid(string gridname,bool exceptionOnUnknown)
{
	if (ptCurrentGrid != NULL && gridname == ptCurrentGrid->sName)
		return;
	ptCurrentGrid = mgNamedGrids[gridname];
	if (ptCurrentGrid == NULL) {
		if (exceptionOnUnknown)
			throw new ValueException(gridname,"unknown GridBox-object",8099);
		ptCurrentGrid = new GridBox();
		ptCurrentGrid->scnScene = this;
		ptCurrentGrid->sName = gridname;
		mgNamedGrids[gridname] = ptCurrentGrid; 
	}
}

void Scene::selectOutputList(string olname,bool exceptionOnUnknown)
{
	if (ptCurrentOutputList != NULL && olname == ptCurrentOutputList->sName)
		return;
	ptCurrentOutputList = moNamedOutputLists[olname];
	if (ptCurrentOutputList == NULL) {
		if (exceptionOnUnknown)
			throw new ValueException(olname,"unknown OutputList-object",8099);
		ptCurrentOutputList = new OutputList();
		ptCurrentOutputList->ptScene = this;
		ptCurrentOutputList->sName = olname;
		moNamedOutputLists[olname] = ptCurrentOutputList; 
	}
}

void Scene::selectOutputFile(string ofname,bool exceptionOnUnknown)
{
	if (ptCurrentOutputFile != NULL && ofname == ptCurrentOutputFile->sName)
		return;
	ptCurrentOutputFile = moNamedOutputFiles[ofname];
	if (ptCurrentOutputFile == NULL) {
		if (exceptionOnUnknown)
			throw new ValueException(ofname,"unknown OutputFile-object",8099);
		ptCurrentOutputFile = new OutputFile();
		ptCurrentOutputFile->ptScene = this;
		ptCurrentOutputFile->sName = ofname;
		moNamedOutputFiles[ofname] = ptCurrentOutputFile; 
	}
}
