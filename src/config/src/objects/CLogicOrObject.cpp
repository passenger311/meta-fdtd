#include "CLogicOrObject.h"

CLogicOrObject::CLogicOrObject(CObject* obj1, CObject* obj2)
{
  oObject1 = obj1;	
  oObject2 = obj2;	
}

CObject* CLogicOrObject::clone() {
  return new CLogicOrObject(oObject1->clone(), oObject2->clone());
}

int CLogicOrObject::Inside(pointframe& pfrm)
{
  if (!box.frameIntersection(pfrm))
    return 0;
  int insd1 = oObject1->Inside(pfrm);
  int insd2 = oObject2->Inside(pfrm);
  return MAX(insd1, insd2);
}
	
bool CLogicOrObject::PointInside(vec3& point)
{
  return oObject1->PointInside(point) || oObject2->PointInside(point);
}
