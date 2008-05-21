#ifndef CCYLINDER_H_
#define CCYLINDER_H_

#include "CObject.h"

/*
 * 
 * 
 */
class CCylinder : public CObject
{
public:
	CCylinder()
	{
		m_dRadius = 1;
		m_dRadiusSquare = 1;
		m_dHalfHeight = 0.5;
	}
	CCylinder(double radius, double height)
	{
		m_dRadius = radius;
		m_dRadiusSquare = radius*radius;
		m_dHalfHeight = height / 2;
	}
	virtual ~CCylinder()
	{ }
	
	virtual void preProcess()
	{
		vec3 p_start(-m_dRadius,-m_dRadius,-m_dHalfHeight); 
		vec3 p_end(m_dRadius,m_dRadius,m_dHalfHeight); 
		box = frame(p_start,p_end);
	}

	virtual bool PointInside(vec3& point);

protected:
	double m_dRadius;
	double m_dRadiusSquare;
	double m_dHalfHeight;
};

#endif /*CCYLINDER_H_*/
