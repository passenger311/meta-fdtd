#include "CEllipsoid.h"

CEllipsoid::CEllipsoid(vec3& pos, double sx, double sy, double sz)
{
  vPosition = pos;
  dScaleX = sx;
  dScaleY = sy;
  dScaleZ = sz;
}
	
CObject* CEllipsoid::clone() {
  return new CEllipsoid(vPosition,dScaleX,dScaleY,dScaleZ);
}

void CEllipsoid::preProcess()
{
  vec3 spos(vPosition-vec3(dScaleX, dScaleY, dScaleZ));
  box = frame(spos, dScaleX*2, dScaleY*2, dScaleZ*2);
}
	
bool CEllipsoid::PointInside(vec3& point)
{
  vec3 dst((point[VX]-vPosition[VX])/dScaleX,
	   (point[VY]-vPosition[VY])/dScaleY,
	   (point[VZ]-vPosition[VZ])/dScaleZ);
  return (dst.length2() <= 1);
}
