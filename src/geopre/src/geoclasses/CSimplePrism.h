#ifndef CSIMPLEPRISM_H_
#define CSIMPLEPRISM_H_

#include "CObject.h"
#include "CLogicAndNotObject.h"
#include "CLogicOrContainer.h"
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
	CSimplePrism()
	{

	}
	CSimplePrism(vector<vec2>& pointlist, double height=1)
	{
		m_vPointlist.clear();
		for (vector<vec2>::iterator pt = pointlist.begin(); pt != pointlist.end(); pt++)
			m_vPointlist.push_back(*pt);
		m_dHalfHeight = height / 2;
	}

	virtual ~CSimplePrism()
	{
		
	}
	
	virtual void preProcess()
	{
		//TODO: Error System
		//TODO: Checking for Intersections and convexity
		//TODO: Simple contained Objects for performance (inner spheres & boxes) 
		if (m_vPointlist.size() < 3)
			throw int(775); //Not enough Points
		vector<vec2> convexShape;
		vector<vec2> tempShape;
		vector<CObject*> subObjects;
		int iMinX = 0;
		for (unsigned int i = 1; i < m_vPointlist.size(); i++)
			if (m_vPointlist[i][VX] < m_vPointlist[iMinX][VX])
				iMinX = i;
		convexShape.push_back(m_vPointlist[iMinX]);
		for (unsigned int i = 1; i < m_vPointlist.size(); i++)
		{
			unsigned int I = (i+iMinX)%m_vPointlist.size();
			vec2 curPt = m_vPointlist[I];
			line2 tmpLine = line2(convexShape.back(), curPt);
			int tmpSgn = 0;
			bool skipPoint = false;
			for (unsigned int k = 1; k < m_vPointlist.size()-MAX(tempShape.size(),1); k++) {
				unsigned int tmpI = (k+i+iMinX)%m_vPointlist.size();
				double tmpDist = tmpLine.signedDistance(m_vPointlist[tmpI]);
				int tmp2Sgn = SIGN(tmpDist);
				if (k == 1 || tmpSgn == 0) {
					tmpSgn = tmp2Sgn;
				} else if ((tmp2Sgn == 1 || tmp2Sgn == -1) && tmpSgn != tmp2Sgn) {
						skipPoint = true;
						break;
				}
			}
			if (skipPoint) {
				if (tempShape.size() == 0)
					tempShape.push_back(convexShape.back());
				tempShape.push_back(curPt);	
			} else {
				if (tmpSgn == -1)
					tmpLine.reverse();
				if (tempShape.size() > 0) {
					// Sub Object
					tempShape.push_back(curPt);	
					subObjects.push_back(new CSimplePrism(tempShape, m_dHalfHeight*2)); 
					tempShape.clear();
				}
				convexShape.push_back(curPt);
			}
		}
		CConvexPrism* outerPrism = new CConvexPrism(convexShape, m_dHalfHeight*2);
		if (subObjects.size() == 0)
			oObject = outerPrism;
		else if (subObjects.size() == 1) {
			CLogicAndNotObject* myobj = new CLogicAndNotObject(outerPrism,subObjects[0]); 
			oObject = myobj;
		} else {
			CLogicOrContainer* subCont = new CLogicOrContainer(subObjects);
			CLogicAndNotObject* myobj = new CLogicAndNotObject(outerPrism,subCont); 
			oObject = myobj;
		}
		oObject->preProcess();
		box = oObject->box;
		
	}
	
	virtual int Inside(pointframe& pfrm)
	{
		return oObject->Inside(pfrm);
	}	
	virtual bool PointInside(vec3& point)
	{
		return oObject->PointInside(point);
	}
};

#endif /*CSIMPLEPRISM_H_*/
