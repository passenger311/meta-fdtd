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
	CBezierPrism()
	{

	}
	CBezierPrism(vector<vec2>& pointlist, double height=1, int steps=10)
	{
		for (vector<vec2>::iterator pt = pointlist.begin(); pt != pointlist.end(); pt++)
			m_vBezierPoints.push_back(*pt);
		m_dHalfHeight = height / 2;
		m_dSteps = steps;
	}
	
	virtual void preProcess()
	{
		m_vPointlist.clear();
		if (m_vBezierPoints.size() % 3 == 0)
			m_vBezierPoints.push_back(m_vBezierPoints.back());
		// TODO: Exception Handling
		if (m_vBezierPoints.size() < 4 || m_vBezierPoints.size() % 3 != 1)
			throw int(458); // Invalid Point Count
		for (vector<vec2>::iterator pt = m_vBezierPoints.begin(); (pt+1) != m_vBezierPoints.end(); pt+=3) {
			vec2 p1 = *pt;
			vec2 p2 = *(pt+1);
			vec2 p3 = *(pt+2);
			vec2 p4 = *(pt+3);
			for (int i = 0; i<m_dSteps; i++) {
				double t = ((double)i)/m_dSteps;
				vec2 r = p1 + 3*t*(p2-p1) + 3*t*t*(p1-2*p2+p3) + t*t*t*(3*p2-3*p3-p1+p4);
				m_vPointlist.push_back(r);
			}
		}
		m_vPointlist.push_back(m_vBezierPoints.back());
		CSimplePrism::preProcess();
	}
};

#endif /*CBEZIERPRISM_H_*/
