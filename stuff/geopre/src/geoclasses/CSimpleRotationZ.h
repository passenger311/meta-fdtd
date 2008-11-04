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
	CSimpleRotationZ()
	{
		oObject = NULL;
	}
	CSimpleRotationZ(CObject* object)
	{
		oObject = object;
	}
	virtual ~CSimpleRotationZ()
	{
		
	}

	virtual bool addSubObject(CObject* object)
	{ 
		if (oObject == NULL)
			oObject = object;
		else
			return false;
		return true;
	}
	
	virtual void preProcess()
	{
		oObject->preProcess();
		double maxX = MAX(fabs(oObject->box.position_start[VX]),fabs(oObject->box.position_end[VX]));
		double maxZ = MAX(fabs(oObject->box.position_start[VZ]),fabs(oObject->box.position_end[VZ]));
		box.position_start = vec3(-maxX,-maxX,-maxZ);
		box.position_end = vec3(maxX,maxX,maxZ);
		box.calcDimensions();
	}
	virtual int Inside(pointframe& pfrm)
	{
		if (!box.frameIntersection(pfrm))
			return 0;
		//return 1;
		pointframe pff;
		pff = pfrm;
		for (int i = 0; i < pff.m_iPointCount; i++) {
			pff.m_vPoints[i][VX] = sqrt(pff.m_vPoints[i][VX]*pff.m_vPoints[i][VX]+pff.m_vPoints[i][VY]*pff.m_vPoints[i][VY]);
			pff.m_vPoints[i][VY] = 0;
		}
		pff.readFrameFromPoints();
		return oObject->Inside(pff); 
	}	
	virtual bool PointInside(vec3& point)
	{
		vec3 pt = point;
		pt[VX] = sqrt(point[VX]*point[VX]+point[VY]*point[VY]);
		pt[VY] = 0;
		return oObject->PointInside(pt);
	}
};

#endif /*CSIMPLEROTATIONZ_H_*/
