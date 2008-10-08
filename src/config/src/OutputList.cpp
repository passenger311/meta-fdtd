#include "OutputList.h"

OutputList::OutputList()
{
	sParameter = "";
	sGrid = "";
	sFormat = "";
	sObjectCut = "";
	bYeeGrid_isset = false;
	bYeeGrid = true;
	bUseCluster_isset = false;
	bUseCluster = false;
	
	ptScene = NULL;	 
	sName = "";
}

OutputList::~OutputList()
{
}

void OutputList::generateOutput()
{
	for (vector<OutputFile*>::iterator opiter = vOutputFiles.begin();
			 opiter != vOutputFiles.end(); opiter++) {
		(*opiter)->generateOutput();		 	
	}
}
