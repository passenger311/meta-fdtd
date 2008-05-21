#include "ArgumentReader.h"

ArgumentReader::ArgumentReader()
{
}

ArgumentReader::~ArgumentReader()
{
}

double ArgumentReader::readDouble(const char* strval , double *defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("","no value given", 2004);
		return *defaultval;
	}
	return atof(strval); // TODO: units	
}

int ArgumentReader::readInteger(const char* strval , int *defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2004);
		return *defaultval;
	}
	return atoi(strval);	
}

bool ArgumentReader::readBoolean(const char* strval, bool *defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2004);
		return *defaultval;
	}
	if (strcmp(strval,"TRUE") == 0)
		return true;
	else if (strcmp(strval,"FALSE") == 0)
		return false;
	else
		throw new ValueParseException(strval,"unknown boolean value, use TRUE or FALSE",2006);	
}

string ArgumentReader::readFilename(const char* strval , string *defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2004);
		return string(*defaultval);
	}
	//TODO: check if valid file
	return string(strval);	
}

string ArgumentReader::readIdentifier(const char* strval , string *defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2004);
		return string(*defaultval);
	}
	//TODO: check if valid identifier
	return string(strval);	
}

vec3 ArgumentReader::readVector3(const char* strval, const vec3* defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2004);
		return vec3(*defaultval);
	}
	char delims[] = ",";
	char *result = NULL;
	char vstrcpy[50]; // TODO: Max Vector String length
	memset( vstrcpy, '\0', 50 );
	strcpy(vstrcpy,strval);
	result = strtok( vstrcpy, delims );
	char *cchar = &vstrcpy[0];
	double vn[3];
	int ix = 0;
	while( ix < 3 && result != NULL ) {
		readDouble(result,&vn[ix]);
		cchar += strlen(result)+1;
		vn[ix] = this->readDouble(result);
		result = strtok( cchar, delims );
		ix++;
	}
	if (ix < 3) {
		throw new ValueParseException(strval, "invalid vector3 string, wroing number of components", 2001);
	}
	return vec3(vn[VX],vn[VY],vn[VZ]);
} 

vec2 ArgumentReader::readVector2(const char* strval, const vec2* defaultval)
{
	if (strval == NULL || strval[0] == 0) {
		if (defaultval == NULL)
			throw new ValueException("", "no value given", 2005);
		return vec2(*defaultval);
	}
	char delims[] = ",";
	char *result = NULL;
	char vstrcpy[50]; // TODO: Max Vector String length
	memset( vstrcpy, '\0', 50 );
	strcpy(vstrcpy,strval);
	char *cchar = &vstrcpy[0];
	result = strtok( vstrcpy, delims );
	double vn[2];
	int ix = 0;
	while( ix < 2 && result != NULL ) {
		readDouble(result,&vn[ix]);
		cchar += strlen(result)+1;
		vn[ix] = this->readDouble(result);
		result = strtok( cchar, delims );
		ix++;
	}
	if (ix < 2) {
		throw new ValueParseException(strval, "invalid vector2 string, wroing number of components", 2002);
	}
	return vec2(vn[VX],vn[VY]);
} 

vector<vec2>* ArgumentReader::readPointlist2(const char* vstr)
{
	vector<vec2>* ret = new vector<vec2>();
	//char delims[] = ";";
	char *result = NULL;
	char vstrcpy[2000]; //TODO: remove length constraint
	memset( vstrcpy, '\0', 2000 );
	strcpy(vstrcpy,vstr);
	char *cchar = &vstrcpy[0];
	result = strtok( vstrcpy, ";:#|" );
	while( result != NULL ) {
		vec2 vn = readVector2(result,NULL);
		ret->push_back(vn);
		cchar += strlen(result)+1;
		result = strtok( cchar, ";:#|" );
	}
	return ret;
}

CObject* ArgumentReader::readObjectReference(const char* vstr)
{
	if (ptScene->moNamedObjects[vstr] == NULL) {
		throw new ValueException(vstr, "unknown object reference", 2003	);
	}
	return ptScene->moNamedObjects[vstr];
}
