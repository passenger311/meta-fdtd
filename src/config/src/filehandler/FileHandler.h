#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "../GridBox.h"

using namespace std;

class GridBox;

class FileHandler
{
protected:
	ofstream m_fsFileStream;
public:
	string sFile;
	GridBox* gbGridBox;
	bool bYeeGrid;
	
public:
	FileHandler();
	virtual ~FileHandler();
	virtual void openFile(string file = "");

	virtual void writeFileHeader();

	virtual void writeGridZDataSlice(int z);
//	virtual void writeGridDataSlice(double* grid);

	virtual void writeFileFooter();

};

#endif /*FILEHANDLER_H_*/
