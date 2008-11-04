#include "Scene.h"
#include "xmlparser/SetupFileParser.h"
#include "GridBox.h"
#include "constants.h"
#include "errorclasses.h"
#include <sstream>

Scene::Scene()
{
	sSceneFile = "";
	sCurrentClusterNode = "";
	bSilentMode = false;
	bBackupOldFiles = false;
	bOverwriteFiles = false;
	bNoErrorOutput = false;
	bFullStatistics = false;
	bQuit = false;
	ptReader = new ArgumentReader();
	ptReader->ptScene = this;
	ptCurrentOutputFile = NULL;
	ptCurrentOutputList = NULL;
	ptCurrentGrid = NULL;
	selectGrid("default");
	selectOutputList("default");
	bSingleOutputFile = false;
	//selectOutputFile("default");
	exFileReading = NULL;
}

Scene::~Scene()
{
}

void Scene::readPriorityCommandLineArguments(int argc, char **argv)
{
	if (argc < 2) {
		throw new Exception("mandatory command-line parameters are missing",1003);	
	}
	int cur_arg = 1;
	while (cur_arg < argc)
	{
		if (strcmp(argv[cur_arg],"--help") == 0) {
			//TODO
			cout << "HELP";
			bQuit = true;
			return;
		} else if (strcmp(argv[cur_arg],"--extended-help") == 0) {
			//TODO
			cout << "EXTENDED HELP";
			bQuit = true;
			return;
		} else if (strcmp(argv[cur_arg],"--license") == 0) {
			//TODO
			cout << PRODUCT_NAME << "  " 
				<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION 
				<< " (" << VERSION_DATE << ")\n";
			// output license.txt
			bQuit = true;
			return;
		} else if (strcmp(argv[cur_arg],"--version") == 0) {
			//TODO
			cout << PRODUCT_NAME << "  " 
				<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION 
				<< " (" << VERSION_DATE << ")\n";
			cout << PRODUCT_COPYRIGHT << "\n";
			cout << PRODUCT_INFO << "\n";
			cout << "Author: " << AUTHOR_NAME << " (" << AUTHOR_MAIL << ")\n";
			bQuit = true;
			return;
		} else if (strcmp(argv[cur_arg],"-s") == 0) {
			bSilentMode = true;
		} else if (strcmp(argv[cur_arg],"-b") == 0) {
			bBackupOldFiles = true;
		} else if (strcmp(argv[cur_arg],"-o") == 0) {
			bOverwriteFiles = true;
		} else if (strcmp(argv[cur_arg],"-n") == 0) {
			bNoErrorOutput = true;
		} else if (strcmp(argv[cur_arg],"-f") == 0) {
			bFullStatistics = true;
		} else if (strcmp(argv[cur_arg],"-e") == 0) {
			Exception::bShowErrorDetails = true;
		}
		cur_arg++;
	}
	sSceneFile = ptReader->readFilename(argv[argc-1]); // TODO: Check file for validity
}

void Scene::readCommandLineArguments(int argc, char **argv)
{
	if (argc < 2) {
		throw new Exception("mandatory command-line parameters are missing",1003);	
	}
	int cur_arg = 1;
	selectGrid("default");
	selectOutputFile("default");
	selectOutputList("default");
	while (cur_arg < argc)
	{
		if (strcmp(argv[cur_arg],"-file") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -file missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			if (ptCurrentOutputFile == NULL)
				ptCurrentOutputFile = new OutputFile();
			ptCurrentOutputFile->sFile = ptReader->readFilename(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-parameter") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -parameter missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			ptCurrentOutputFile->sParameter = ptReader->readIdentifier(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-format") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -format missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			ptCurrentOutputFile->sFormat = ptReader->readIdentifier(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-grid") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -grid missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			ptCurrentOutputFile->sGrid = ptReader->readIdentifier(argv[++cur_arg]);
			selectGrid(ptCurrentOutputFile->sGrid); 
		} else if (strcmp(argv[cur_arg],"-objectcut") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -objectcut missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			ptCurrentOutputFile->sObjectCut = ptReader->readIdentifier(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-oyee") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -objectcut missing",1002);	
			bSingleOutputFile = true;
			selectOutputFile("default");
			ptCurrentOutputFile->bYeeGrid = ptReader->readBoolean(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-points") == 0) {
			if (cur_arg + 4 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -points missing",1002);	
			ptCurrentGrid->iCellsX = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iCellsY = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iCellsZ = ptReader->readInteger(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-origin") == 0) {
			if (cur_arg + 4 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -origin missing",1002);
			vec3 vc;	
			vc[VX] = ptReader->readDouble(argv[++cur_arg]); 
			vc[VY] = ptReader->readDouble(argv[++cur_arg]); 
			vc[VZ] = ptReader->readDouble(argv[++cur_arg]);
			ptCurrentGrid->frBBox.position_end = ptCurrentGrid->frBBox.position_end 
				- ptCurrentGrid->frBBox.position_start + vc;
			ptCurrentGrid->frBBox.position_start = vc;
		} else if (strcmp(argv[cur_arg],"-size") == 0) {
			if (cur_arg + 4 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -size missing",1002);
			vec3 vc;	
			vc[VX] = ptReader->readDouble(argv[++cur_arg]); 
			vc[VY] = ptReader->readDouble(argv[++cur_arg]); 
			vc[VZ] = ptReader->readDouble(argv[++cur_arg]);
			ptCurrentGrid->frBBox.position_end = ptCurrentGrid->frBBox.position_start + vc; 
		} else if (strcmp(argv[cur_arg],"-sg") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -sg missing",1002);	
			string subgridding = ptReader->readIdentifier(argv[++cur_arg]);
			if (subgridding == "TRUE") {
				ptCurrentGrid->bNoSubgridding = false;
				ptCurrentGrid->bAlwaysSubgridding = false;
			}else if (subgridding == "FALSE") {
				ptCurrentGrid->bNoSubgridding = true;
				ptCurrentGrid->bAlwaysSubgridding = false;
			}else if (subgridding == "ALWAYS") {
				ptCurrentGrid->bNoSubgridding = false;
				ptCurrentGrid->bAlwaysSubgridding = true;
			} else
				throw new ValueParseException(argv[cur_arg-2],"unknown subgridding value, use TRUE, FALSE or ALWAYS",1004);	
		} else if (strcmp(argv[cur_arg],"-yee") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -yee missing",1002);	
			ptCurrentGrid->bYeeGrid = ptReader->readBoolean(argv[++cur_arg]);
		} else if (strcmp(argv[cur_arg],"-sg_divisions") == 0) {
			if (cur_arg + 4 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -sg_divisions missing",1002);	
			ptCurrentGrid->iSubGriddingDivX = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iSubGriddingDivY = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iSubGriddingDivZ = ptReader->readInteger(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-pf_scale") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -pf_scale missing",1002);	
			ptCurrentGrid->dPointframeScale = ptReader->readDouble(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-pf_divisions") == 0) {
			if (cur_arg + 4 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -pf_divisions missing",1002);	
			ptCurrentGrid->iPointframeDivisionsX = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iPointframeDivisionsY = ptReader->readInteger(argv[++cur_arg]); 
			ptCurrentGrid->iPointframeDivisionsZ = ptReader->readInteger(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-pf_facepoints") == 0) {
			if (cur_arg + 2 >=  argc)
				throw new ArgumentException(argv[cur_arg],"parameter(s) for -pf_facepoints missing",1002);	
			ptCurrentGrid->bPointframeFacepoints = ptReader->readBoolean(argv[++cur_arg]); 
		} else if (strcmp(argv[cur_arg],"-cn") == 0) {
			if (cur_arg + 2 ==  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -cn missing",1002);	
			sCurrentClusterNode = argv[++cur_arg];
		} else if (strcmp(argv[cur_arg],"-om") == 0) {
			if (cur_arg + 2 ==  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -om missing",1002);	
			selectOutputList(argv[++cur_arg],true);
		} else if (strcmp(argv[cur_arg],"-os") == 0) {
			if (cur_arg + 2 ==  argc)
				throw new ArgumentException(argv[cur_arg],"parameter for -os missing",1002);	
			selectOutputFile(argv[++cur_arg],true);
			bSingleOutputFile = true;
		} else if (strcmp(argv[cur_arg],"-s") == 0) {
			bSilentMode = true;
		} else if (strcmp(argv[cur_arg],"-b") == 0) {
			bBackupOldFiles = true;
		} else if (strcmp(argv[cur_arg],"-o") == 0) {
			bOverwriteFiles = true;
		} else if (strcmp(argv[cur_arg],"-n") == 0) {
			bNoErrorOutput = true;
		} else if (strcmp(argv[cur_arg],"-f") == 0) {
			bFullStatistics = true;
		} else if (strcmp(argv[cur_arg],"-e") == 0) {
			Exception::bShowErrorDetails = true;
		} else {
			if (cur_arg + 1 != argc)
				throw new ArgumentException(argv[cur_arg],"unknown parameter",1001);	
		}
		cur_arg++;
	}
	
}

void Scene::readSceneFile()
{
	SetupFileParser parser(this);
	const char *filename;
	FILE* xmlFile;
	char buf[BUFSIZ];
	int done;
	selectGrid("default");
	filename = sSceneFile.c_str();
	if (strlen(filename)==0)
		throw new FileException(filename, "no filename given",3001);
	
	xmlFile = fopen(filename, "r");
	if (!xmlFile)
		throw new FileException(filename, "could not open file",3002);
	
	do {
		size_t len = fread(buf, 1, sizeof(buf), xmlFile /*stdin*/);
		done = len < sizeof(buf);
		if (!parser.XML_Parse(buf, len, done)) {
			stringstream sstr;
			sstr << "XML parsing error: " 
				<< XML_ErrorString(parser.XML_GetErrorCode())
				<< " at line "
				<< parser.XML_GetCurrentLineNumber();
			throw new FileParsingException(filename, sstr.str(), 3003);
		}
		if (exFileReading != NULL) {
			throw exFileReading;
		}
	} while (!done);
}

void Scene::buildScene()
{
	
}

void Scene::generateOutput()
{
	//selectOutputFile("default",true);
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
