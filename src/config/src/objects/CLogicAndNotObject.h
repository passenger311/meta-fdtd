#ifndef CLOGICANDNOTOBJECT_H_
#define CLOGICANDNOTOBJECT_H_

#include "CBinaryObject.h"

class CLogicAndNotObject : public CBinaryObject
{
public:
  CLogicAndNotObject(CObject* obj1, CObject* obj2);
  virtual CObject* clone();
  virtual int Inside(pointframe& pfrm);
  virtual void preProcess();
  virtual bool PointInside(vec3& point);
};

#endif /*CLOGICANDNOTOBJECT_H_*/
