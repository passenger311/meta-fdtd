#include "CLogicAndNotObject.h"

CLogicAndNotObject::CLogicAndNotObject(CObject* obj1, CObject* obj2)
{
  oObject1 = obj1;	
  oObject2 = obj2;	
}

CObject* CLogicAndNotObject::clone() {
  return new CLogicAndNotObject(oObject1->clone(), oObject2->clone());
}

int CLogicAndNotObject::Inside(pointframe& pfrm)
{
  if (!box.frameIntersection(pfrm))
    return 0;
  int insd1 = oObject1->Inside(pfrm);
  int insd2 = oObject2->Inside(pfrm);
  if (insd1 == 2 && insd2 == 0)
    return 2;
  if ((insd1 == 1 && (insd2 == 0 || insd2 == 1)) ||
      (insd1 == 2 && insd2 == 1))
    return 1;
  return 0;
}

void CLogicAndNotObject::preProcess()
{
  oObject1->preProcess();
  oObject2->preProcess();
  box = oObject1->box;
}

bool CLogicAndNotObject::PointInside(vec3& point)
{
  return oObject1->PointInside(point) && !oObject2->PointInside(point);
}
