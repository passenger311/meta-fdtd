#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "../Grid.h"

using namespace std;

class Grid;

class FileHandler
{
protected:
  ofstream m_fsFileStream;
public:
  string sFile;

public:
  FileHandler();
  FileHandler(const char* fname);
  virtual ~FileHandler();
  virtual void writeFileHeader(Grid* gbGrid);
  virtual void writeGridZDataSlice(Grid* gbGrid, int z);
  virtual void writeFileFooter(Grid* gbGrid);
};

#endif /*FILEHANDLER_H_*/
