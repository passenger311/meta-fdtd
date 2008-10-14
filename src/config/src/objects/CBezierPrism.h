#ifndef CBEZIERPRISM_H_
#define CBEZIERPRISM_H_

#include "CSimplePrism.h"

class CBezierPrism : public CSimplePrism
{
protected:
	vector<vec2> m_vBezierPoints;
	int m_dSteps;

public:

public:
	CBezierPrism();
	CBezierPrism(vector<vec2>& pointlist, double height=1, int steps=10);
	virtual CObject* clone();
	virtual void preProcess();
};

#endif /*CBEZIERPRISM_H_*/
