#ifndef GRID_H_
#define GRID_H_

#include "objects/algebra3extension.h"
#include "filehandler/FileHandler.h"
#include "Scene.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

enum YEEGRIDCOMPONENT {ZERO, EX, EY, EZ, HX, HY, HZ}; //YeeGrid components

class FileHandler;
class Scene;

class Grid
{
public:
	frame frBBox;
	int iCellsX, iCellsY, iCellsZ;
	double dPointframeScale;
	int iPointframeDivisionsX,iPointframeDivisionsY,iPointframeDivisionsZ;
	bool bPointframeFacepoints;
	int iSubGriddingDivX, iSubGriddingDivY, iSubGriddingDivZ;
	bool bAlwaysSubgridding;
	bool bNoSubgridding;
	bool bYeeGrid;
	Scene* scnScene;
	int iSubGriddedCells,
		iSubGriddedFilledCells,
		iSubGriddedEmptyCells,
		iFilledGridCells,
		iEmptyGridCells;
	string sName;

protected:
	int m_iCellsRealX, m_iCellsRealY, m_iCellsRealZ;
	double *m_dGridvalues;
	double *m_dGridvalues_1;
	double *m_dGridvalues_2;
	int m_iCurrentZ;
public:
	Grid();
	Grid(Grid* gridbox);
	Grid(Scene* scene, frame box, int cellsX, int cellsY, int cellsZ);

	~Grid();
	
	void generateOutput(FileHandler* fhd);
	
	double getDataPoint(int x, int y, int z, YEEGRIDCOMPONENT cp);
};

#endif /*GRID_H_*/
