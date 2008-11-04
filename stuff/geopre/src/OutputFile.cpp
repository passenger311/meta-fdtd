#include "OutputFile.h"
#include "errorclasses.h"
#include "filehandler/FileHandlerVTK.h"
#include "filehandler/FileHandlerFortranIN.h"

OutputFile::OutputFile()
{
	sFile = "";
	sParameter = "";
	sGrid = "default";
	sFormat = "";
	sObjectCut = "";
	bYeeGrid = true;
	bUseCluster = false;
}

OutputFile::~OutputFile()
{
}

void OutputFile::generateOutput()
{
	ptScene->selectGrid(sGrid,true);
	FileHandler* fhd;
	if (sFormat == "" || sFormat == "default" || sFormat == "vtk") {
		fhd = new FileHandlerVTK();
		sFormat = "vtk";
	} else if (sFormat == "fortran_in")
		fhd = new FileHandlerFortranIN();
	else
		throw new ValueParseException(sFormat,"unknown file format",7012);
	fhd->sFile = sFile;
	fhd->bYeeGrid = bYeeGrid;
	fhd->gbGridBox = ptScene->ptCurrentGrid;
	
	// TODO: parameter selection & cluster crop
	if (!ptScene->bSilentMode) {
		cout << "Generating file \"" << sFile << "\"" << endl; 	
		cout << "-> Format: " << sFormat << endl; 	
		cout << "-> Grid: " << sGrid << endl; 	
		cout << "-> Parameter: " << sParameter << endl; 	
		cout << "generating.";
	}
	GridBox* cg = ptScene->ptCurrentGrid;
	cg->generateOutput(fhd);
	if (!ptScene->bSilentMode) {
		cout << "...finished!" << endl;
		cout << "* Subgridded Cells: " << cg->iSubGriddedCells << "\n";
		cout << "* Subgridded filled Cells: " << cg->iSubGriddedFilledCells << "\n";
		cout << "* Subgridded empty Cells: " << cg->iSubGriddedEmptyCells << "\n";
		cout << "* Subgridded partly filled Cells: " << (cg->iSubGriddedCells-cg->iSubGriddedEmptyCells-cg->iSubGriddedFilledCells) << "\n";
		cout << "* Filled Cells: " << cg->iFilledGridCells << "\n";
		cout << "* Empty Cells: " << cg->iEmptyGridCells << "\n";
		cout << endl;
	}

	delete fhd;
}
