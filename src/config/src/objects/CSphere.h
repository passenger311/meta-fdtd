#ifndef CSPHERE_H_
#define CSPHERE_H_

#include "CObject.h"

class CSphere : public CObject
{

protected:
    double m_dRadiusSqr;

public:
    vec3 vPosition;
    double dRadius;

public:
	
    CSphere(vec3& pos, double radius);
    virtual CObject* clone();
    virtual void preProcess();
    virtual bool PointInside(vec3& point);
};

#endif /*CSPHERE_H_*/
