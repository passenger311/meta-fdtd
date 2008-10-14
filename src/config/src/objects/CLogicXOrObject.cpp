#include "CLogicXOrObject.h"

CLogicXOrObject::CLogicXOrObject(CObject* obj1, CObject* obj2)
{
  oObject1 = obj1;	
  oObject2 = obj2;	
}

CObject* CLogicXOrObject::clone() {
  return new CLogicXOrObject(oObject1->clone(),oObject2->clone());
}

int CLogicXOrObject::Inside(pointframe& pfrm)
{
  if (!box.frameIntersection(pfrm))
    return 0;
  int insd1 = oObject1->Inside(pfrm);
  int insd2 = oObject2->Inside(pfrm);
  if (insd1 == 0 && insd2 == 0)
    return 0;
  if ((insd1 == 0 && insd2 == 2) ||
      (insd1 == 2 && insd2 == 0))
    return 2;
  return 1;
}

bool CLogicXOrObject::PointInside(vec3& point)
{
  return oObject1->PointInside(point) ^ oObject2->PointInside(point);
}
