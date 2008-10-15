
#include "CSphere.h"


CSphere::CSphere(vec3& pos, double radius) {
  vPosition = pos;
  dRadius = radius;
#ifdef DEBUG
  printf("Sphere: at = %f %f %f radius = %f\n",vPosition[0],vPosition[1],vPosition[2],dRadius);
#endif
}

CObject* CSphere::clone()
{
  return new CSphere(vPosition, dRadius);
}


void CSphere::preProcess()
{
  vec3 spos(vPosition-vec3(dRadius, dRadius, dRadius));
  box = frame(spos, dRadius*2, dRadius*2, dRadius*2);
  m_dRadiusSqr = dRadius*dRadius;	
}
	
bool CSphere::PointInside(vec3& point)
{
  vec3 dst(point[VX]-vPosition[VX],point[VY]-vPosition[VY],point[VZ]-vPosition[VZ]);
  return (dst.length2() <= m_dRadiusSqr);
}


