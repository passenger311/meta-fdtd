#include "FileHandler.h"

FileHandler::FileHandler()
{
}

FileHandler::FileHandler(const char* fname)
{
  sFile = fname;
}

FileHandler::~FileHandler()
{
}

void FileHandler::writeFileHeader(Grid* gbGrid)
{
  m_fsFileStream.open(sFile.c_str());
}

void FileHandler::writeGridZDataSlice(Grid* gbGrid, int z)
{

}

void FileHandler::writeFileFooter(Grid* gbGrid)
{
  m_fsFileStream.close();
}
