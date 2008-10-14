#ifndef CBOX_H_
#define CBOX_H_

#include "CObject.h"

class CBox : public CObject
{
 public:
	
  CBox(vec3& pos, vec3& size);
  CBox(frame& ff);
  virtual CObject* clone();	
  virtual bool PointInside(vec3& point);
};

#endif /*CBOX_H_*/
