#ifndef ERRORCLASSES_H_
#define ERRORCLASSES_H_

#include <string>
#include <sstream>

using namespace std;

class Exception
{
public:
	Exception();
	Exception(string msg, int errnum);
		
	virtual ~Exception();
	
	virtual string getMessage();

public:
	static bool bShowErrorDetails;
	
protected:
	string m_sMessage;
	int m_iErrorNumber;
	
};

class ArgumentException : public Exception
{
public:
	ArgumentException();
	ArgumentException(string argumentname, string msg, int errnum);
	
	virtual ~ArgumentException() {  }
	
	virtual string getMessage();
	
protected:
	string m_sArgument;
};

class ValueParseException : public Exception
{
public:
	ValueParseException();
	ValueParseException(string argumentname, string msg, int errnum);
	
	virtual ~ValueParseException() {  }
	
	virtual string getMessage();
	
protected:
	string m_sArgument;
};

class ValueException : public Exception
{
public:
	ValueException();
	ValueException(string argumentname, string msg, int errnum);
	
	virtual ~ValueException() {  }
	
	virtual string getMessage();
	
protected:
	string m_sArgument;
};

class FileException : public Exception
{
public:
	FileException();
	FileException(string argumentname, string msg, int errnum);
	
	virtual ~FileException() {  }
	
	virtual string getMessage();
	
protected:
	string m_sArgument;
};

class FileParsingException : public Exception
{
public:
	FileParsingException();
	FileParsingException(string argumentname, string msg, int errnum);
	
	virtual ~FileParsingException() {  }
	
	virtual string getMessage();
	
protected:
	string m_sArgument;
};


#endif /*ERRORCLASSES_H_*/
