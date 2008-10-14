#ifndef CLOGICXOROBJECT_H_
#define CLOGICXOROBJECT_H_

#include "CBinaryObject.h"

class CLogicXOrObject : public CBinaryObject
{
public:

  CLogicXOrObject(CObject* obj1, CObject* obj2);
  virtual CObject* clone();
  virtual int Inside(pointframe& pfrm);
  virtual bool PointInside(vec3& point);
};

#endif /*CLOGICXOROBJECT_H_*/
