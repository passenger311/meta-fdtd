#ifndef FILEHANDLERVTK_H_
#define FILEHANDLERVTK_H_

#include "FileHandler.h"

class FileHandlerVTK : public FileHandler
{
public:
  FileHandlerVTK();
  FileHandlerVTK(const char* fname);
  virtual ~FileHandlerVTK();
  virtual void writeFileHeader(Grid* gbGrid);
  virtual void writeGridZDataSlice(Grid* gbGrid, int z);
};

#endif /*FILEHANDLERVTK_H_*/
