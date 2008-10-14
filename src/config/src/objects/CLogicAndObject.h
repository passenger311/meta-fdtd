#ifndef CLOGICANDOBJECT_H_
#define CLOGICANDOBJECT_H_

#include "CBinaryObject.h"

class CLogicAndObject : public CBinaryObject
{
public:
  CLogicAndObject(CObject* obj1, CObject* obj2);
  CObject* clone();
  virtual int Inside(pointframe& pfrm);
  virtual void preProcess();
  virtual bool PointInside(vec3& point);
};

#endif /*CLOGICANDOBJECT_H_*/
