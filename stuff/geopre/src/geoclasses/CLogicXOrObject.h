#ifndef CLOGICXOROBJECT_H_
#define CLOGICXOROBJECT_H_

class CLogicXOrObject : public CBinaryObject
{
public:
	CLogicXOrObject(CObject* obj1, CObject* obj2)
	{
		oObject1 = obj1;	
		oObject2 = obj2;	
	}
	virtual ~CLogicXOrObject()
	{ }
	virtual int Inside(pointframe& pfrm)
	{
		if (!box.frameIntersection(pfrm))
			return 0;
		int insd1 = oObject1->Inside(pfrm);
		int insd2 = oObject2->Inside(pfrm);
		if (insd1 == 0 && insd2 == 0)
			return 0;
		if ((insd1 == 0 && insd2 == 2) ||
			(insd1 == 2 && insd2 == 0))
			return 2;
		return 1;
	}
	
	virtual bool PointInside(vec3& point)
	{
		return oObject1->PointInside(point) ^ oObject2->PointInside(point);
	}
};

#endif /*CLOGICXOROBJECT_H_*/
