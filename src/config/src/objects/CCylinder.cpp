#include "CCylinder.h"


CCylinder::CCylinder()
{
  dRadius = 1;
  m_dRadiusSquare = 1;
  m_dHalfHeight = 0.5;
}
CCylinder::CCylinder(vec3& pos, double radius, double height)
{
  vPosition = pos;
  dRadius = radius;
  m_dRadiusSquare = radius*radius;
  m_dHalfHeight = height / 2.0;
}

CObject* CCylinder::clone() {
  return new CCylinder(vPosition, dRadius, m_dHalfHeight*2.0);
}
	
void CCylinder::preProcess()
{
  vec3 v(dRadius,dRadius,m_dHalfHeight);
  vec3 p_start(vPosition - v); 
  vec3 p_end(vPosition + v); 
  box = frame(p_start,p_end);
}

bool CCylinder::PointInside(vec3& point)
{
  vec3 p = point - vPosition;
  return p[VX]*p[VX]+p[VY]*p[VY] <= m_dRadiusSquare &&
    -m_dHalfHeight <= p[VZ] && m_dHalfHeight >= p[VZ];
}
