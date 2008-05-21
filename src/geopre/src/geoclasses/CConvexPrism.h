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
	CConvexPrism()
	{

	}
	CConvexPrism(vector<vec2>& pointlist, double height=1)
	{
		m_vPointlist.clear();
		for (vector<vec2>::iterator pt = pointlist.begin(); pt != pointlist.end(); pt++)
			m_vPointlist.push_back(*pt);
		m_dHalfHeight = height / 2;
	}

	virtual ~CConvexPrism()
	{
		
	}
	
	virtual void preProcess()
	{
		//TODO: Error System
		//TODO: Checking for Intersections and convexity
		//TODO: Simple contained Objects for performance (inner spheres & boxes) 
		if (m_vPointlist.size() < 3)
			throw int(775); //Not enough Points
		vec3* pl = new vec3[m_vPointlist.size()];
		int i = 0;
		for (vector<vec2>::iterator pt = m_vPointlist.begin(); pt != m_vPointlist.end(); pt++)
		{
			pl[i] = (*pt);
			if (i != 0) {
				line2 ln(*(pt-1),*(pt));
				if (i == 1 && ln.signedDistance(*(pt+1)) < 0)
					ln.reverse();
				if (i > 1 && ln.signedDistance(*(pt-2)) < 0)
					ln.reverse();
				m_lBorders.push_back(ln);		
			} else {
				line2 ln(m_vPointlist.back(),*pt);
				if (ln.signedDistance(*(pt+1)) < 0)
					ln.reverse();
				m_lBorders.push_back(ln);
			}
			i++;
		}
		if (i != (int) m_vPointlist.size())
			throw int(775);
		box.getFromPointlist(pl,i);
		box.position_start[VZ] = -m_dHalfHeight;
		box.position_end[VZ] = m_dHalfHeight;
		box.calcDimensions();
		delete pl;
		
	}
	
	//TODO: critical Points and dimensions
//	virtual int Inside(pointframe& pfrm)
//	{
//		return 2; 
//	}	
	virtual bool PointInside(vec3& point)
	{
		for (vector<line2>::iterator iline = m_lBorders.begin(); iline != m_lBorders.end(); iline++)
		{
			vec2 v2(point,VZ);
			if (iline->signedDistance(v2) < 0)
				return false;
		}
		return true;
	}
};

#endif /*CCONVEXPRISM_H_*/
