#ifndef CLOGICOROBJECT_H_
#define CLOGICOROBJECT_H_

#include "CBinaryObject.h"

class CLogicOrObject : public CBinaryObject
{
public:
  CLogicOrObject(CObject* obj1, CObject* obj2);
  virtual CObject * clone();
  virtual int Inside(pointframe& pfrm);
  virtual bool PointInside(vec3& point);
};

#endif /*CLOGICOROBJECT_H_*/
