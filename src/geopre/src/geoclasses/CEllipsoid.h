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
//	CBox();
	
	virtual ~CEllipsoid()
	{ }
	
	CEllipsoid(vec3& pos, double sx, double sy, double sz)
	{
		vPosition = pos;
		dScaleX = sx;
		dScaleY = sy;
		dScaleZ = sz;
	}
	
	virtual void preProcess()
	{
		vec3 spos(vPosition-vec3(dScaleX, dScaleY, dScaleZ));
		box = frame(spos, dScaleX*2, dScaleY*2, dScaleZ*2);
	}
	
	virtual bool PointInside(vec3& point)
	{
		vec3 dst((point[VX]-vPosition[VX])/dScaleX,
			(point[VY]-vPosition[VY])/dScaleY,
			(point[VZ]-vPosition[VZ])/dScaleZ);
		return (dst.length2() <= 1);
	}
};

#endif /*CELLIPSOID_H_*/
