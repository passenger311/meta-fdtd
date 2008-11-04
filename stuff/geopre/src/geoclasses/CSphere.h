#ifndef CSPHERE_H_
#define CSPHERE_H_

#include "CObject.h"

class CSphere : public CObject
{
protected:
	double m_dRadiusSqr;
public:
	vec3 vPosition;
	double dRadius;
public:
//	CBox();
	
	virtual ~CSphere()
	{ }
	
	CSphere(vec3& pos, double radius)
	{
		vPosition = pos;
		dRadius = radius;
	}
	
	virtual void preProcess()
	{
		vec3 spos(vPosition-vec3(dRadius, dRadius, dRadius));
		box = frame(spos, dRadius*2, dRadius*2, dRadius*2);
		m_dRadiusSqr = dRadius*dRadius;	
	}
	
	virtual bool PointInside(vec3& point)
	{
		vec3 dst(point[VX]-vPosition[VX],point[VY]-vPosition[VY],point[VZ]-vPosition[VZ]);
		return (dst.length2() <= m_dRadiusSqr);
	}
};

#endif /*CSPHERE_H_*/
