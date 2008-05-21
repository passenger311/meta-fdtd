#ifndef ARGUMENTREADER_H_
#define ARGUMENTREADER_H_

#include <string>
#include <vector>
#include "geoclasses/algebra3extension.h"
#include "geoclasses/CObject.h"
#include "errorclasses.h"
#include "Scene.h"

using namespace std;

class Scene;

class ArgumentReader
{

public:
	ArgumentReader();
	virtual ~ArgumentReader();
	
	double readDouble(const char* strval, double *defaultval = NULL);
	int readInteger(const char* strval, int *defaultval = NULL);
	bool readBoolean(const char* strval, bool *defaultval = NULL);
	string readFilename(const char* strval, string *defaultval = NULL);
	string readIdentifier(const char* strval , string *defaultval = NULL);
	vec3 readVector3(const char* strval, const vec3 *defaultval = NULL);
	vec2 readVector2(const char* strval, const vec2 *defaultval = NULL);
	vector<vec2>* readPointlist2(const char* vstr);
	CObject* readObjectReference(const char* vstr);

	Scene* ptScene;

};

#endif /*ARGUMENTREADER_H_*/
