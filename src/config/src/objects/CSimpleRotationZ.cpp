#include "CSimpleRotationZ.h"

CSimpleRotationZ::CSimpleRotationZ()
{
  oObject = NULL;
}

CSimpleRotationZ::CSimpleRotationZ(CObject* object)
{
  oObject = object;
}

CObject* CSimpleRotationZ::clone() {
  return new CSimpleRotationZ(oObject->clone());
}

bool CSimpleRotationZ::addSubObject(CObject* object)
{ 
  if (oObject == NULL)
    oObject = object;
  else
    return false;
  return true;
}
	
void CSimpleRotationZ::preProcess()
{
  oObject->preProcess();
  double maxX = MAX(fabs(oObject->box.position_start[VX]),fabs(oObject->box.position_end[VX]));
  double maxZ = MAX(fabs(oObject->box.position_start[VZ]),fabs(oObject->box.position_end[VZ]));
  box.position_start = vec3(-maxX,-maxX,-maxZ);
  box.position_end = vec3(maxX,maxX,maxZ);
  box.calcDimensions();
}

int CSimpleRotationZ::Inside(pointframe& pfrm)
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

bool CSimpleRotationZ::PointInside(vec3& point)
{
  vec3 pt = point;
  pt[VX] = sqrt(point[VX]*point[VX]+point[VY]*point[VY]);
  pt[VY] = 0;
  return oObject->PointInside(pt);
}
