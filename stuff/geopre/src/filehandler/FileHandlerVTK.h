#ifndef FILEHANDLERVTK_H_
#define FILEHANDLERVTK_H_

#include "FileHandler.h"

class FileHandlerVTK : public FileHandler
{
public:
	FileHandlerVTK();
	virtual ~FileHandlerVTK();

	virtual void writeFileHeader();

	virtual void writeGridZDataSlice(int z);
//	virtual void writeGridDataSlice(double* grid);
};

#endif /*FILEHANDLERVTK_H_*/
