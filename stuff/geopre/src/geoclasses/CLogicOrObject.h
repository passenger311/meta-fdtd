#ifndef CLOGICOROBJECT_H_
#define CLOGICOROBJECT_H_

#include "CBinaryObject.h"

class CLogicOrObject : public CBinaryObject
{
public:
	CLogicOrObject(CObject* obj1, CObject* obj2)
	{
		oObject1 = obj1;	
		oObject2 = obj2;	
	}
	virtual ~CLogicOrObject();
	virtual int Inside(pointframe& pfrm)
	{
		if (!box.frameIntersection(pfrm))
			return 0;
		int insd1 = oObject1->Inside(pfrm);
		int insd2 = oObject2->Inside(pfrm);
		return MAX(insd1, insd2);
	}
	
	virtual bool PointInside(vec3& point)
	{
		return oObject1->PointInside(point) || oObject2->PointInside(point);
	}
};

#endif /*CLOGICOROBJECT_H_*/
