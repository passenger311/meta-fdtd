#include "Grid.h"

void Grid::Init() {
  frBBox.position_start = vec3(-1,-1,-1);
  frBBox.position_end = vec3(1,1,1);
  iCellsX = 10;
  iCellsY = 10;
  iCellsZ = 10;
  dPointframeScale = 2;
  iPointframeDivisionsX = 1;
  iPointframeDivisionsY = 1;
  iPointframeDivisionsZ = 1;
  bPointframeFacepoints = true;
  iSubGriddingDivX = 5;
  iSubGriddingDivY = 5;
  iSubGriddingDivZ = 5;
  bAlwaysSubgridding = false;
  bNoSubgridding = false;
  bYeeGrid = true;
  m_iCurrentZ = -1;
}

Grid::Grid()
{
  Init();
}

Grid::Grid(Grid* grid) {
  Init();
  frBBox = grid->frBBox;
  iCellsX = grid->iCellsX;
  iCellsY = grid->iCellsY;
  iCellsZ = grid->iCellsZ;
  dPointframeScale = grid->dPointframeScale;
  iPointframeDivisionsX = grid->iPointframeDivisionsX;
  iPointframeDivisionsY = grid->iPointframeDivisionsY;
  iPointframeDivisionsZ = grid->iPointframeDivisionsZ;
  bPointframeFacepoints = grid->bPointframeFacepoints;
  iSubGriddingDivX = grid->iSubGriddingDivX;
  iSubGriddingDivY = grid->iSubGriddingDivY;
  iSubGriddingDivZ = grid->iSubGriddingDivZ;
  bAlwaysSubgridding = grid->bAlwaysSubgridding;
  bNoSubgridding = grid->bNoSubgridding;
  bYeeGrid = grid->bYeeGrid;
}

Grid::Grid(frame& box, int cellsX, int cellsY, int cellsZ, int offsX, int offsY, int offsZ)
{
  Init();
  frBBox = box;
  
  iCellsX = cellsX;
  iCellsY = cellsY;
  iCellsZ = cellsZ;
  iOffsetX = offsX;
  iOffsetY = offsY;
  iOffsetZ = offsZ;

  if (cellsX < 0 ) { // interprete frame coordinates as grid coordinates
    double from = frBBox.position_start[VX];
    double to = frBBox.position_end[VX];    
    iCellsX = (int)(to-from+1.5);
    iOffsetX = (int)from;
  } 
  if (cellsY < 0 ) { // interprete frame coordinates as grid coordinates
    double from = frBBox.position_start[VY];
    double to = frBBox.position_end[VY];    
    iCellsY = (int)(to-from+1.5);
    iOffsetY = (int)from;
  } 
  if (cellsZ < 0 ) { // interprete frame coordinates as grid coordinates
    double from = frBBox.position_start[VZ];
    double to = frBBox.position_end[VZ];    
    iCellsZ = (int)(to-from+1.5);
    iOffsetZ = (int)from;
  } 
}

Grid::~Grid()
{
}

void Grid::generateOutput(Scene* scnScene, FileHandler* fhd, const char* method, bool silentMode)
{
  fhd->writeFileHeader(this);
	
  vector<CObject*> objects(scnScene->objects);
  // Reset Statistics
  iSubGriddedCells = 0;
  iSubGriddedFilledCells = 0;
  iSubGriddedEmptyCells = 0;
  iFilledGridCells = 0;
  iEmptyGridCells = 0;
  int cellsX = iCellsX;
  int cellsY = iCellsY;
  int cellsZ = iCellsZ;
  double dDimX = frBBox.position_end[VX]-frBBox.position_start[VX]; 
  double dDimY = frBBox.position_end[VY]-frBBox.position_start[VY];
  double dDimZ = frBBox.position_end[VZ]-frBBox.position_start[VZ];
  vec3 vSMcell = frBBox.position_start;
  if (bYeeGrid) {
    cellsX = iCellsX*2+1;
    cellsY = iCellsY*2+1;
    cellsZ = iCellsZ*2+1;
    if (iCellsX > 1) 
      dDimX = dDimX / (iCellsX - 1) / 2 * cellsX;
    else
      dDimX *= 1.5;
    if (iCellsY > 1) 
      dDimY = dDimY / (iCellsY - 1) / 2 * cellsY;
    else
      dDimY *= 1.5;
    if (iCellsZ > 1) 
      dDimZ = dDimZ / (iCellsZ - 1) / 2 * cellsZ;
    else
      dDimZ *= 1.5;
  } else {
    // Total Dimensions
    if (cellsX > 1) 
      dDimX = dDimX / (cellsX - 1) * cellsX;
    if (cellsY > 1) 
      dDimY = dDimY / (cellsY - 1) * cellsY;
    if (cellsZ > 1) 
      dDimZ = dDimZ / (cellsZ - 1) * cellsZ;
  }
  m_iCellsRealX = cellsX;
  m_iCellsRealY = cellsY;
  m_iCellsRealZ = cellsZ;
  double dDX = dDimX / cellsX;
  double dDY = dDimY / cellsY;
  double dDZ = dDimZ / cellsZ;
  if (bYeeGrid) vSMcell -= 0.5*vec3(dDX,dDY,dDZ);
  // Subgridding values
  double dSubValSteps = 1.0 / (iSubGriddingDivX*iSubGriddingDivY*iSubGriddingDivZ);
  double dSubDX = dDX / iSubGriddingDivX;
  double dSubDY = dDY / iSubGriddingDivY;
  double dSubDZ = dDZ / iSubGriddingDivZ;
  // Datagrid
  m_dGridvalues = new double[cellsX*cellsY];
  m_dGridvalues_1 = bYeeGrid ? new double[cellsX*cellsY] : NULL;
  m_dGridvalues_2 = bYeeGrid ? new double[cellsX*cellsY] : NULL;
  // Setup Pointframe
  vec3 startp(vSMcell - dPointframeScale/2*vec3(dDX,dDY,dDZ));
  vec3 endp(vSMcell + dPointframeScale/2*vec3(dDX,dDY,dDZ));
  pointframe pfrm;
  if (bNoSubgridding)
    pfrm = pointframe(startp,endp,0,0,0);
  else {
    pfrm = pointframe(startp,endp,iPointframeDivisionsX,iPointframeDivisionsY,iPointframeDivisionsZ,bPointframeFacepoints);
    pfrm.movingPoint = vSMcell + vec3(-0.5*dDX,0.5*dDY,0.5*dDZ) + vec3(dSubDX/2,dSubDY/2,dSubDZ/2);	
  } 
  pfrm.moveX(dDimX);
  pfrm.moveY(dDimY);
  for (vector<CObject*>::iterator scnobj = objects.begin(); 
       scnobj != objects.end(); scnobj++)
    {
      (*scnobj)->preProcess();
    }
  // --> [AH]
  // reverse list, so that later defined objects paint over previous defined ones
  std::reverse(objects.begin(), objects.end()); 
  // ... and sort according to depth
  std::sort(objects.begin(), objects.end(), PredicateCObjectDepth());
  // <-- [AH]
  vector<CObject*> currentobjects;

  for (int iZ = 0; iZ < cellsZ; iZ++) {
    int iGrid = 0;
    pfrm.moveY(-dDimY);
    currentobjects.clear();
    for (vector<CObject*>::iterator scnobj = objects.begin(); 
	 scnobj != objects.end(); scnobj++) {
      if (pfrm.frameIntersectionZ((*scnobj)->box)) currentobjects.push_back(*scnobj);
    }
    //		if (currentobjects.size() == 0)
    //			continue;
    for (int iY = 0; iY < cellsY; iY++) {
      pfrm.moveX(-dDimX);
      for (int iX = 0; iX < cellsX; iX++) {

	if (!bNoSubgridding) {

	  // -------------------------

	  double value = -1; // set bad value
	  int insd = 0;
	  for (vector<CObject*>::iterator scnobj = currentobjects.begin(); 
	       scnobj != currentobjects.end(); scnobj++) {
	    insd = MAX(insd, (*scnobj)->iLastInside = (*scnobj)->Inside(pfrm));
	    value = (*scnobj)->dValue; // +[AH]
	    if (insd == 2 && !bAlwaysSubgridding) break;
	  }

	  if (insd == 1 || bAlwaysSubgridding) {

	    iSubGriddedCells++;	
	    vec3 ptSub(pfrm.movingPoint);
	    //ptSub += vec3(iX*dDX+dSubDX/2,(iY+1)*dDY+dSubDY/2,(iZ+1)*dDZ+dSubDZ/2);
	    
	    double dGridVal = 0;
	    for (int iSubX = 0; iSubX < iSubGriddingDivX; iSubX++) {
	      ptSub[VY] -= dDY; 
	      for (int iSubY = 0; iSubY < iSubGriddingDivY; iSubY++) {
		ptSub[VZ] -= dDZ; 
		for (int iSubZ = 0; iSubZ < iSubGriddingDivZ; iSubZ++) {
		  bool subinside = false;
		  for (vector<CObject*>::iterator scnobj = currentobjects.begin(); 
		       scnobj != currentobjects.end(); scnobj++) {
		    if ((*scnobj)->iLastInside > 0 && (*scnobj)->PointInside(ptSub)) {
		      subinside = true;
		      value = (*scnobj)->dValue;
		      break;
		    }
		  }
		  dGridVal += subinside ? value : scnScene->dValue ;
		  ptSub[VZ] += dSubDZ; 
		}
		ptSub[VY] += dSubDY; 
	      }
	      ptSub[VX] += dSubDX; 
	    }
	    dGridVal *= dSubValSteps;
	    m_dGridvalues[iGrid] = dGridVal;
	    if (dGridVal  >= 0.99999999) iSubGriddedFilledCells++;	
	    if (dGridVal == 0.0) iSubGriddedEmptyCells++;	
	
	  } else {

	    m_dGridvalues[iGrid] = (insd == 2) ? value : scnScene->dValue;
	    
	  }

	} else {

	  // ---- no subgridding -----

 	  m_dGridvalues[iGrid] = scnScene->dValue; // background value
	  for (vector<CObject*>::iterator scnobj = currentobjects.begin(); 
	       scnobj != currentobjects.end(); scnobj++) {
	    if ((*scnobj)->InsideFrame(pfrm.movingPoint) &&  
		(*scnobj)->PointInside(pfrm.movingPoint)) {
	      m_dGridvalues[iGrid] = (*scnobj)->dValue; // overwrite background
	      break;
	    }
	  }

	  // -------------------------

	}


	// ---> [AH]
	// FIX: the gridvalues can be > 1 if the value is > 1.
	if (m_dGridvalues[iGrid] >= 0.99999999 && m_dGridvalues[iGrid] <= 1.0000001 ) {
	  m_dGridvalues[iGrid] = 1;
	}
	
	if (m_dGridvalues[iGrid] >= 0.99999999) { iFilledGridCells++; }
	
	// if (m_dGridvalues[iGrid] >= 0.99999999) {
	//	m_dGridvalues[iGrid] = 1;
	//	iFilledGridCells++;	
	// }
	// <---
	
	if (m_dGridvalues[iGrid] == 0.0) iEmptyGridCells++;
	iGrid++;
	pfrm.moveX(dDX);
      }
      pfrm.moveY(dDY);
    }
    if (bYeeGrid) {
      if (iZ >= 2 && iZ % 2 == 0) {
	m_iCurrentZ = iZ / 2 - 1;
	fhd->writeGridZDataSlice(this,m_iCurrentZ);
	// fhd->writeGridDataSlice(m_dGridvalues);
      }
      double *gvtemp = m_dGridvalues_2;
      m_dGridvalues_2 = m_dGridvalues_1;
      m_dGridvalues_1 = m_dGridvalues;
      m_dGridvalues = gvtemp;
    } else {
      m_iCurrentZ = iZ;
      fhd->writeGridZDataSlice(this,m_iCurrentZ);
    }
    if (!silentMode && iZ % MAX((int)cellsZ/50,1) == 0) {
      cout << "." << flush;
    }
    pfrm.moveZ(dDZ);
  }
  delete m_dGridvalues;
  if (bYeeGrid) {
    delete m_dGridvalues_1;
    delete m_dGridvalues_2;
  }
  fhd->writeFileFooter(this);
}


double Grid::getDataPoint(int x, int y, int z, YEEGRIDCOMPONENT cp)
{
  if (z != m_iCurrentZ)
    throw int(8734); //TODO: Error Handling
  if (bYeeGrid) {
    int baseindex = (2*y+1)*m_iCellsRealX + (2*x+1);
    switch (cp) {
    case EX:
      return (m_dGridvalues_2[baseindex - m_iCellsRealX] + // (2x+1,2y,2z)
	      m_dGridvalues_2[baseindex - m_iCellsRealX+1] + // (2x+2,2y,2z)
	      m_dGridvalues_2[baseindex] + // (2x+1,2y+1,2z)
	      m_dGridvalues_2[baseindex + 1] + // (2x+2,2y+1,2z)
	      m_dGridvalues_1[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+1)
	      m_dGridvalues_1[baseindex - m_iCellsRealX + 1] + // (2x+2,2y,2z+1)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + 1])/8; // (2x+2,2y+1,2z+1)
    case EY:
      return (m_dGridvalues_2[baseindex - 1] + // (2x,2y+1,2z)
	      m_dGridvalues_2[baseindex] + // (2x+1,2y+1,2z)
	      m_dGridvalues_2[baseindex + m_iCellsRealX - 1] + // (2x,2y+2,2z)
	      m_dGridvalues_2[baseindex + m_iCellsRealX] + // (2x+1,2y+2,2z)
	      m_dGridvalues_1[baseindex - 1] + // (2x,2y+1,2z+1)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX - 1] + // (2x,2y+2,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX])/8; // (2x+1,2y+2,2z+1)
    case EZ:
      return (m_dGridvalues_1[baseindex - m_iCellsRealX - 1] + // (2x,2y,2z+1)
	      m_dGridvalues_1[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+1)
	      m_dGridvalues_1[baseindex - 1] + // (2x,2y+1,2z+1)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues[baseindex - m_iCellsRealX - 1] + // (2x,2y,2z+2)
	      m_dGridvalues[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+2)
	      m_dGridvalues[baseindex - 1] + // (2x,2y+1,2z+2)
	      m_dGridvalues[baseindex])/8; // (2x+1,2y+1,2z+2)
    case HX:
      return (m_dGridvalues_1[baseindex - 1] + // (2x,2y+1,2z+1)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX - 1] + // (2x,2y+2,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX] + // (2x+1,2y+2,2z+1)
	      m_dGridvalues[baseindex - 1] + // (2x,2y+1,2z+2)
	      m_dGridvalues[baseindex] + // (2x+1,2y+1,2z+2)
	      m_dGridvalues[baseindex + m_iCellsRealX - 1] + // (2x,2y+2,2z+2)
	      m_dGridvalues[baseindex + m_iCellsRealX])/8; // (2x+1,2y+2,2z+2)
    case HY:
      return (m_dGridvalues_1[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+1)
	      m_dGridvalues_1[baseindex - m_iCellsRealX + 1] + // (2x+2,2y,2z+1)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + 1] + // (2x+2,2y+1,2z+1)
	      m_dGridvalues[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+2)
	      m_dGridvalues[baseindex - m_iCellsRealX + 1] + // (2x+2,2y,2z+2)
	      m_dGridvalues[baseindex] + // (2x+1,2y+1,2z+2)
	      m_dGridvalues[baseindex + 1])/8; // (2x+2,2y+1,2z+2)
    case HZ:
      return (m_dGridvalues_2[baseindex] + // (2x+1,2y+1,2z)
	      m_dGridvalues_2[baseindex + 1] + // (2x+2,2y+1,2z)
	      m_dGridvalues_2[baseindex + m_iCellsRealX] + // (2x+1,2y+2,2z)
	      m_dGridvalues_2[baseindex + m_iCellsRealX + 1] + // (2x+2,2y+2,2z)
	      m_dGridvalues_1[baseindex] + // (2x+1,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + 1] + // (2x+2,2y+1,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX] + // (2x+1,2y+2,2z+1)
	      m_dGridvalues_1[baseindex + m_iCellsRealX + 1])/8; // (2x+2,2y+2,2z+1)
    default:
      return (m_dGridvalues_2[baseindex - m_iCellsRealX - 1] + // (2x,2y,2z)
	      m_dGridvalues_2[baseindex - m_iCellsRealX] + // (2x+1,2y,2z)
	      m_dGridvalues_2[baseindex - 1] + // (2x,2y+1,2z)
	      m_dGridvalues_2[baseindex] + // (2x+1,2y+1,2z)
	      m_dGridvalues_1[baseindex - m_iCellsRealX - 1] + // (2x,2y,2z+1)
	      m_dGridvalues_1[baseindex - m_iCellsRealX] + // (2x+1,2y,2z+1)
	      m_dGridvalues_1[baseindex - 1] + // (2x,2y+1,2z+1)
	      m_dGridvalues_1[baseindex])/8; // (2x+1,2y+1,2z+1)
    }
  } else
    return m_dGridvalues[y*m_iCellsRealX + x];
}
