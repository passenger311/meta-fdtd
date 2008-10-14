#ifndef CBINARYOBJECT_H_
#define CBINARYOBJECT_H_

#include "CObject.h"

class CBinaryObject : public CObject
{
public:
	CObject* oObject1;
	CObject* oObject2;
	
public:
	CBinaryObject();
	CBinaryObject(CObject* obj1, CObject* obj2);
	virtual void preProcess()
	{
		oObject1->preProcess();
		oObject2->preProcess();
		box = oObject1->box.combineWith(oObject2->box);
	}
	virtual int Inside(pointframe& pfrm)
	{
		return 0;
	}
	
	virtual bool PointInside(vec3& point)
	{
		return false;
	}
	virtual bool addSubObject(CObject* object)
	{ 
		if (oObject1 == NULL)
			oObject1 = object;
		else if (oObject2 == NULL)
			oObject2 = object;
		else
			return false;
		return true;
	}

};

#endif /*CBINARYOBJECT_H_*/
