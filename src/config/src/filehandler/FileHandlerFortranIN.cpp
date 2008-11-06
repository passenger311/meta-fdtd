#include "FileHandlerFortranIN.h"

FileHandlerFortranIN::FileHandlerFortranIN()
{
}

FileHandlerFortranIN::FileHandlerFortranIN(const char* fname, int comps) : FileHandler(fname)
{
  m_iComps = comps;
}

FileHandlerFortranIN::~FileHandlerFortranIN()
{
}

void FileHandlerFortranIN::writeFileHeader(Grid* gbGrid)
{
  double from[3];
  from[0] = gbGrid->frBBox.position_start[VX];
  from[1] = gbGrid->frBBox.position_start[VY];
  from[2] = gbGrid->frBBox.position_start[VZ];
  double to[3];
  to[0] = gbGrid->frBBox.position_end[VX];
  to[1] = gbGrid->frBBox.position_end[VY];
  to[2] = gbGrid->frBBox.position_end[VZ];
  double 
    sizeX = to[0]-from[0], 
    sizeY = to[1]-from[1],
    sizeZ = to[2]-from[2];
  int 
    i0 = gbGrid->iOffsetX,
    j0 = gbGrid->iOffsetY,
    k0 = gbGrid->iOffsetZ;
  int
    i1 = i0 + gbGrid->iCellsX -1,
    j1 = j0 + gbGrid->iCellsY -1,
    k1 = k0 + gbGrid->iCellsZ -1;


  double dDX = sizeX / (gbGrid->iCellsX-1);
  if (gbGrid->iCellsX <= 1) dDX = sizeX;
  double dDY = sizeY / (gbGrid->iCellsY-1);
  if (gbGrid->iCellsY <= 1) dDY = sizeY;
  double dDZ = sizeZ / (gbGrid->iCellsZ-1);
  if (gbGrid->iCellsZ <= 1) dDZ = sizeZ;
  int iPointCount = gbGrid->iCellsX * gbGrid->iCellsY * gbGrid->iCellsZ;
  
  m_fsFileStream.open(sFile.c_str());
  
  m_fsFileStream << " ! (.IN) FIELD DATA FILE" << endl;
  m_fsFileStream << " ! NUMNODES: " << iPointCount << endl;
  m_fsFileStream << "(SET" << endl;
  m_fsFileStream << i0 << " " << i1 << " " << 1 << " "
		 << j0 << " " << j1 << " " << 1 << " "
		 << k0 << " " << k1 << " " << 1 << endl;
}

void FileHandlerFortranIN::writeGridZDataSlice(Grid* gbGrid,int z)
{
  double gridData;
  for (int y=0; y<gbGrid->iCellsY; y++) {
    for (int x=0; x<gbGrid->iCellsX; x++) {
      if (gbGrid->bYeeGrid) {
	m_fsFileStream
	  << "  "
	  << gbGrid->getDataPoint(x,y,z,EX) << ' '
	  << gbGrid->getDataPoint(x,y,z,EY) << ' '
	  << gbGrid->getDataPoint(x,y,z,EZ);
	if ( m_iComps == 3 ) { 
	  m_fsFileStream << "\n"; 
	} else {
	  m_fsFileStream << ' ' 
	  << gbGrid->getDataPoint(x,y,z,HX) << ' '
	  << gbGrid->getDataPoint(x,y,z,HY) << ' '
	  << gbGrid->getDataPoint(x,y,z,HZ) << "\n";
	}
      } else {
        gridData = gbGrid->getDataPoint(x,y,z,ZERO); 
	m_fsFileStream
	  << "  "
	  << gridData << ' '
	  << gridData << ' '
	  << gridData;
	if ( m_iComps == 3 ) { 
	  m_fsFileStream << "\n"; 
	} else {
	  m_fsFileStream << ' '
	  << gridData << ' '
	  << gridData << ' '
	  << gridData << "\n";
	}
      }
    }
  }
}

void FileHandlerFortranIN::writeFileFooter(Grid* gbGrid)
{
  m_fsFileStream << ")SET" << endl;
  m_fsFileStream.close();
}
