#ifndef CLOGICANDNOTOBJECT_H_
#define CLOGICANDNOTOBJECT_H_

#include "CBinaryObject.h"

class CLogicAndNotObject : public CBinaryObject
{
public:
	CLogicAndNotObject(CObject* obj1, CObject* obj2)
	{
		oObject1 = obj1;	
		oObject2 = obj2;	
	}
	virtual ~CLogicAndNotObject()
	{ }
	virtual int Inside(pointframe& pfrm)
	{
		if (!box.frameIntersection(pfrm))
			return 0;
		int insd1 = oObject1->Inside(pfrm);
		int insd2 = oObject2->Inside(pfrm);
		if (insd1 == 2 && insd2 == 0)
			return 2;
		if ((insd1 == 1 && (insd2 == 0 || insd2 == 1)) ||
			(insd1 == 2 && insd2 == 1))
			return 1;
		return 0;
	}

	virtual void preProcess()
	{
		oObject1->preProcess();
		oObject2->preProcess();
		box = oObject1->box;
	}
	
	virtual bool PointInside(vec3& point)
	{
		return oObject1->PointInside(point) && !oObject2->PointInside(point);
	}
};

#endif /*CLOGICANDNOTOBJECT_H_*/
