#include "CObject.h"

CObject::CObject()
{
	name = "";
	fWeight = 1.;
	iDepth = 0;
}

CObject::~CObject()
{
}

void CObject::preProcess()
{
	
}
int CObject::Inside(pointframe& pfrm)
{	
	if (!box.frameIntersection(pfrm))
		return 0;
	int iInsidePoints = 0;
	for (int i = 0; i < pfrm.m_iPointCount; i++) {
		if (PointInside(pfrm.m_vPoints[i]))
		{
			iInsidePoints++;
			if (iInsidePoints <= i)
				return 1;
		} 
	}
	return (iInsidePoints == pfrm.m_iPointCount) ? 2 : 0;	
}

inline bool CObject::PointInside(vec3& point)
{
	return false;	
}

inline bool CObject::InsideFrame(vec3& point)
{
	return box.pointInFrame(point);	
}

inline bool CObject::InsideFrame(pointframe& pfrm)
{
	return pfrm.frameIntersection(box);
}
