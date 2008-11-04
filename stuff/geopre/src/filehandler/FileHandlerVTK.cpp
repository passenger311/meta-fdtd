#include "FileHandlerVTK.h"

FileHandlerVTK::FileHandlerVTK()
{
}

FileHandlerVTK::~FileHandlerVTK()
{
}

void FileHandlerVTK::writeFileHeader()
{
	double dimX = gbGridBox->frBBox.position_end[VX]-gbGridBox->frBBox.position_start[VX], 
		dimY = gbGridBox->frBBox.position_end[VY]-gbGridBox->frBBox.position_start[VY], 
		dimZ = gbGridBox->frBBox.position_end[VZ]-gbGridBox->frBBox.position_start[VZ];
	double dDX = dimX / (gbGridBox->iCellsX-1);
	if (gbGridBox->iCellsX <= 1)
		dDX = dimX;
	double dDY = dimY / (gbGridBox->iCellsY-1);
	if (gbGridBox->iCellsY <= 1)
		dDY = dimY;
	double dDZ = dimZ / (gbGridBox->iCellsZ-1);
	if (gbGridBox->iCellsZ <= 1)
		dDZ = dimZ;
	int iPointCount = gbGridBox->iCellsX * gbGridBox->iCellsY * gbGridBox->iCellsZ;
	m_fsFileStream << "# vtk DataFile Version 2.0\n";
	m_fsFileStream << "GridSetup generated file\n";
	m_fsFileStream << "ASCII\n";
	m_fsFileStream << "DATASET STRUCTURED_POINTS\n";
	m_fsFileStream << "DIMENSIONS " << gbGridBox->iCellsX << " " << gbGridBox->iCellsY << " " << gbGridBox->iCellsZ << "\n";
	m_fsFileStream << "SPACING " << dDX << " " << dDY << " " << dDZ << "\n";
	m_fsFileStream << "ORIGIN " << gbGridBox->frBBox.position_start[VX] 
			<< " " << gbGridBox->frBBox.position_start[VY]
			<< " " << gbGridBox->frBBox.position_start[VZ] << "\n";
	m_fsFileStream << "\n";
	m_fsFileStream << "POINT_DATA " << iPointCount << "\n";
	if (bYeeGrid)
		m_fsFileStream << "SCALARS values double 6\n";
	else
		m_fsFileStream << "SCALARS values double 1\n";
	m_fsFileStream << "LOOKUP_TABLE default\n";
}

void FileHandlerVTK::writeGridZDataSlice(int z)
{
	for (int y=0; y<gbGridBox->iCellsY; y++) {
		for (int x=0; x<gbGridBox->iCellsX; x++) {
			if (bYeeGrid)
				m_fsFileStream
//						 << x << ':' << y << ':' << z << ' '
						 << gbGridBox->getDataPoint(x,y,z,EX) << ' '
						 << gbGridBox->getDataPoint(x,y,z,EY) << ' '
						 << gbGridBox->getDataPoint(x,y,z,EZ) << ' '
						 << gbGridBox->getDataPoint(x,y,z,HX) << ' '
						 << gbGridBox->getDataPoint(x,y,z,HY) << ' '
						 << gbGridBox->getDataPoint(x,y,z,HZ) << "\n";
			else
				m_fsFileStream 
//					<< x << ':' << y << ':' << z << ' '
					<< gbGridBox->getDataPoint(x,y,z,ZERO) << "\n";
		}
	}
}
