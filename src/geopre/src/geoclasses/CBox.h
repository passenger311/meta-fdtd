#ifndef CBOX_H_
#define CBOX_H_

#include "CObject.h"

class CBox : public CObject
{
public:
//	CBox();
	
	virtual ~CBox()
	{ }
	
	CBox(frame& ff)
	{
		box = ff;
	}
	
	virtual bool PointInside(vec3& point)
	{
		return box.pointInFrame(point);
	}
};

#endif /*CBOX_H_*/
