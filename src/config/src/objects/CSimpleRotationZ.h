#ifndef CSIMPLEROTATIONZ_H_
#define CSIMPLEROTATIONZ_H_

#include "CObject.h"
#include <cmath>

using namespace std;

class CSimpleRotationZ : public CObject
{
protected:

public:

	CObject* oObject;
	
 public:
	CSimpleRotationZ();
	CSimpleRotationZ(CObject* object);
	virtual CObject* clone();
	virtual bool addSubObject(CObject* object);
	virtual void preProcess();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
};

#endif /*CSIMPLEROTATIONZ_H_*/
