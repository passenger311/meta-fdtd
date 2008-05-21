#ifndef COBJECT_H_
#define COBJECT_H_

#include "algebra3extension.h"
#include <string>

using namespace std;

class CObject
{
public:
	CObject();
	virtual ~CObject();
	

public:
	frame box;
	int iLastInside;
	string name;
	
	virtual void preProcess();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
	virtual bool InsideFrame(vec3& point);
	virtual bool InsideFrame(pointframe& pfrm);
	virtual bool addSubObject(CObject* object)
	{ return false; }
};

#endif /*COBJECT_H_*/
