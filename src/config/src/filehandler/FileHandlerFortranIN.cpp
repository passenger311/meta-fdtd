#include "FileHandlerFortranIN.h"

FileHandlerFortranIN::FileHandlerFortranIN()
{
}

FileHandlerFortranIN::~FileHandlerFortranIN()
{
}

void FileHandlerFortranIN::writeFileHeader()
{
	double dimX = gbGrid->frBBox.position_end[VX]-gbGrid->frBBox.position_start[VX], 
		dimY = gbGrid->frBBox.position_end[VY]-gbGrid->frBBox.position_start[VY], 
		dimZ = gbGrid->frBBox.position_end[VZ]-gbGrid->frBBox.position_start[VZ];
	double dDX = dimX / (gbGrid->iCellsX-1);
	if (gbGrid->iCellsX <= 1)
		dDX = dimX;
	double dDY = dimY / (gbGrid->iCellsY-1);
	if (gbGrid->iCellsY <= 1)
		dDY = dimY;
	double dDZ = dimZ / (gbGrid->iCellsZ-1);
	if (gbGrid->iCellsZ <= 1)
		dDZ = dimZ;
	int iPointCount = gbGrid->iCellsX * gbGrid->iCellsY * gbGrid->iCellsZ;
	m_fsFileStream << " ! (SET-) FIELD DATA FILE" << endl;
	m_fsFileStream << " ! META: xxx" << endl;
	m_fsFileStream << " ! ISBOX: xxx" << endl;
	m_fsFileStream << " ! NUMNODES: " << iPointCount << endl;
	m_fsFileStream << " ! IRANGE: xxx" << endl;
	m_fsFileStream << " ! JRANGE: xxx" << endl;
	m_fsFileStream << " ! KRANGE: xxx" << endl;
	m_fsFileStream << "(SET" << endl;
//	m_fsFileStream << "ASCII\n";
//	m_fsFileStream << "DATASET STRUCTURED_POINTS\n";
//	m_fsFileStream << "DIMENSIONS " << gbGrid->iCellsX << " " << gbGrid->iCellsY << " " << gbGrid->iCellsZ << "\n";
//	m_fsFileStream << "SPACING " << dDX << " " << dDY << " " << dDZ << "\n";
//	m_fsFileStream << "ORIGIN " << gbGrid->frBBox.position_start[VX] 
//			<< " " << gbGrid->frBBox.position_start[VY]
//			<< " " << gbGrid->frBBox.position_start[VZ] << "\n";
//	m_fsFileStream << "\n";
//	m_fsFileStream << "POINT_DATA " << iPointCount << "\n";
//	if (bYeeGrid)
//		m_fsFileStream << "SCALARS values double 6\n";
//	else
//		m_fsFileStream << "SCALARS values double 1\n";
//	m_fsFileStream << "LOOKUP_TABLE default\n";
}

void FileHandlerFortranIN::writeGridZDataSlice(int z)
{
	for (int y=0; y<gbGrid->iCellsY; y++) {
		for (int x=0; x<gbGrid->iCellsX; x++) {
			if (bYeeGrid)
				m_fsFileStream
//						 << x << ':' << y << ':' << z << ' '
						 << "  "
						 << gbGrid->getDataPoint(x,y,z,EX) << ' '
						 << gbGrid->getDataPoint(x,y,z,EY) << ' '
						 << gbGrid->getDataPoint(x,y,z,EZ) << ' '
						 << gbGrid->getDataPoint(x,y,z,HX) << ' '
						 << gbGrid->getDataPoint(x,y,z,HY) << ' '
						 << gbGrid->getDataPoint(x,y,z,HZ) << "\n";
			else
				m_fsFileStream 
//					<< x << ':' << y << ':' << z << ' '
						 << "  "
						 << gbGrid->getDataPoint(x,y,z,ZERO) << ' '
						 << gbGrid->getDataPoint(x,y,z,ZERO) << ' '
						 << gbGrid->getDataPoint(x,y,z,ZERO) << ' '
						 << gbGrid->getDataPoint(x,y,z,ZERO) << ' '
						 << gbGrid->getDataPoint(x,y,z,ZERO) << ' '
						 << gbGrid->getDataPoint(x,y,z,ZERO) << "\n";
		}
	}
}

void FileHandlerFortranIN::writeFileFooter()
{
	m_fsFileStream << ")SET" << endl;
	m_fsFileStream.close();
}
