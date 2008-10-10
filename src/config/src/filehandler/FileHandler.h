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
	Grid* gbGrid;
	bool bYeeGrid;
	
public:
	FileHandler();
	virtual ~FileHandler();
	virtual void openFile(string file = "");

	virtual void writeFileHeader();

	virtual void writeGridZDataSlice(int z);

	virtual void writeFileFooter();

};

#endif /*FILEHANDLER_H_*/
