#ifndef CCYLINDER_H_
#define CCYLINDER_H_

#include "CObject.h"

/*
 * 
 * 
 */
class CCylinder : public CObject
{

protected:
  double m_dRadiusSquare;
  double m_dHalfHeight;

public:
  double dRadius;
  vec3 vPosition;

public:
  CCylinder();
  CCylinder(vec3& pos, double radius, double height);
  virtual CObject* clone();
  virtual void preProcess();
  virtual bool PointInside(vec3& point);

};

#endif /*CCYLINDER_H_*/
