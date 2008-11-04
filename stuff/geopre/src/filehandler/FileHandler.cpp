#include "FileHandler.h"

FileHandler::FileHandler()
{
	bYeeGrid = true;
}

FileHandler::~FileHandler()
{
}

void FileHandler::openFile(string file)
{
	if (file != "")
			sFile = file;
	m_fsFileStream.open(sFile.c_str());
}

void FileHandler::writeFileHeader()
{

}

void FileHandler::writeGridZDataSlice(int z)
{

}

void FileHandler::writeFileFooter()
{
	m_fsFileStream.close();
}
