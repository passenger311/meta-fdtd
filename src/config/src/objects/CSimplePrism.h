#ifndef CSIMPLEPRISM_H_
#define CSIMPLEPRISM_H_

#include "CObject.h"
#include "CLogicAndNotObject.h"
#include "CLogicOrContainer.h"
#include "CConvexPrism.h"
#include <vector>

using namespace std;

class CSimplePrism : public CObject
{
protected:
	vector<vec2> m_vPointlist;
	CObject* oObject;
	double m_dHalfHeight;

public:

public:
	CSimplePrism();
	CSimplePrism(vector<vec2>& pointlist, double height=1.);
	virtual CObject* clone();
	virtual void preProcess();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
};

#endif /*CSIMPLEPRISM_H_*/
