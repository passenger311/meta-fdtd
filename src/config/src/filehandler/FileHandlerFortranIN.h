#ifndef FILEHANDLERFORTRANIN_H_
#define FILEHANDLERFORTRANIN_H_

#include "FileHandler.h"

class FileHandlerFortranIN : public FileHandler
{
public:
  FileHandlerFortranIN();
  FileHandlerFortranIN(const char* fname);
  virtual ~FileHandlerFortranIN();
  virtual void writeFileHeader(Grid* gbGrid);
  virtual void writeGridZDataSlice(Grid* gbGrid, int z);
  virtual void writeFileFooter(Grid* gbGrid);
};

#endif /*FILEHANDLERFORTRANIN_H_*/
