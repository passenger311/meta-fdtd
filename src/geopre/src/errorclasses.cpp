
#include "errorclasses.h"

bool Exception::bShowErrorDetails = false;

Exception::Exception()
{
	m_sMessage = "";
}
Exception::Exception(string msg, int errnum = 0)
{
	m_sMessage = msg;	
	m_iErrorNumber = errnum;
}

Exception::~Exception() { }

string Exception::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}


ArgumentException::ArgumentException() 
	: Exception()
{
	m_sArgument = "";
}
ArgumentException::ArgumentException(string argumentname, string msg, int errnum = 0)
	: Exception(msg,errnum)
{
	m_sArgument = argumentname;
}


string ArgumentException::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << " (" << m_sArgument << ")" << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}

ValueException::ValueException() 
	: Exception()
{
	m_sArgument = "";
}
ValueException::ValueException(string argumentname, string msg, int errnum = 0)
	: Exception(msg,errnum)
{
	m_sArgument = argumentname;
}


string ValueException::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << " (" << m_sArgument << ")" << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}

FileException::FileException() 
	: Exception()
{
	m_sArgument = "";
}
FileException::FileException(string argumentname, string msg, int errnum = 0)
	: Exception(msg,errnum)
{
	m_sArgument = argumentname;
}


string FileException::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << " (" << m_sArgument << ")" << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}

FileParsingException::FileParsingException() 
	: Exception()
{
	m_sArgument = "";
}
FileParsingException::FileParsingException(string argumentname, string msg, int errnum = 0)
	: Exception(msg,errnum)
{
	m_sArgument = argumentname;
}


string FileParsingException::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << " (" << m_sArgument << ")" << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}

ValueParseException::ValueParseException() 
	: Exception()
{
	m_sArgument = "";
}
ValueParseException::ValueParseException(string argumentname, string msg, int errnum = 0)
	: Exception(msg,errnum)
{
	m_sArgument = argumentname;
}


string ValueParseException::getMessage()
{
	stringstream tmp;
	tmp << m_sMessage << " (" << m_sArgument << ")" << endl;
	if (Exception::bShowErrorDetails && m_iErrorNumber) {
		tmp << "Error Code: " << m_iErrorNumber << endl;
	}
	return tmp.str();	
}
