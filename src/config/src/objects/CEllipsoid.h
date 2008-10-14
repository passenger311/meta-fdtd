#ifndef CELLIPSOID_H_
#define CELLIPSOID_H_

#include "CObject.h"

class CEllipsoid : public CObject
{
protected:
	double m_dRadiusSqr;
public:
	vec3 vPosition;
	double dScaleX;
	double dScaleY;
	double dScaleZ;
public:
	
	CEllipsoid(vec3& pos, double sx, double sy, double sz);
	virtual CObject* clone();
	virtual void preProcess();
	virtual bool PointInside(vec3& point);
};

#endif /*CELLIPSOID_H_*/
