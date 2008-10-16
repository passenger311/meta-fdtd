#ifndef FILEHANDLERFORTRANIN_H_
#define FILEHANDLERFORTRANIN_H_

#include "FileHandler.h"

class FileHandlerFortranIN : public FileHandler
{
private:
  int m_iComps;
public:
  FileHandlerFortranIN();
  FileHandlerFortranIN(const char* fname, int comps = 6 );
  virtual ~FileHandlerFortranIN();
  virtual void writeFileHeader(Grid* gbGrid);
  virtual void writeGridZDataSlice(Grid* gbGrid, int z);
  virtual void writeFileFooter(Grid* gbGrid);
};

#endif /*FILEHANDLERFORTRANIN_H_*/
