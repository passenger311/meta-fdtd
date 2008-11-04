#ifndef FILEHANDLERFORTRANIN_H_
#define FILEHANDLERFORTRANIN_H_

#include "FileHandler.h"

class FileHandlerFortranIN : public FileHandler
{
public:
	FileHandlerFortranIN();
	virtual ~FileHandlerFortranIN();

	virtual void writeFileHeader();

	virtual void writeGridZDataSlice(int z);
//	virtual void writeGridDataSlice(double* grid);

	virtual void writeFileFooter();
};

#endif /*FILEHANDLERFORTRANIN_H_*/
