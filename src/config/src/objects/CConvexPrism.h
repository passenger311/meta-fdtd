#ifndef CCONVEXPRISM_H_
#define CCONVEXPRISM_H_

#include "CObject.h"
#include <vector>

using namespace std;

class CConvexPrism : public CObject
{
protected:
	vector<vec2> m_vPointlist;
	vector<line2> m_lBorders;
	double m_dHalfHeight;

public:

public:
	CConvexPrism();
	CConvexPrism(vector<vec2>& pointlist, double height=1);
	virtual CObject* clone();
	virtual void preProcess();
//TODO: critical Points and dimensions
//	int Inside(pointframe& pfrm)
//	{
//		return 2; 
//	}	
	virtual bool PointInside(vec3& point);
};

#endif /*CCONVEXPRISM_H_*/
