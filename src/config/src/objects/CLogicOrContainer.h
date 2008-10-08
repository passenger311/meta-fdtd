#ifndef CLOGICORCONTAINER_H_
#define CLOGICORCONTAINER_H_
#include "CObject.h"

class CLogicOrContainer : public CObject
{
protected:
	vector<CObject*> m_vObjectList;	

public:
	CLogicOrContainer(vector<CObject*> objectlist)
	{
		for (vector<CObject*>::iterator iter = objectlist.begin(); iter != objectlist.end(); iter++)
			m_vObjectList.push_back(*iter);
	}
	virtual ~CLogicOrContainer() {  }
	virtual int Inside(pointframe& pfrm)
	{
		if (!box.frameIntersection(pfrm))
			return 0;
		int insd = 0;
		for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++)
			insd = MAX(insd, (*iter)->Inside(pfrm));
		return insd;
	}
	
	virtual bool PointInside(vec3& point)
	{
		for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++)
			if ((*iter)->PointInside(point))
				return true;
		return false;
	}
	virtual bool addSubObject(CObject* object)
	{ 
		m_vObjectList.push_back(object);
		return true;
	}
	virtual void preProcess()
	{
		for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++) {
			(*iter)->preProcess();
			if (iter != m_vObjectList.begin())
				box = box.combineWith((*iter)->box);
			else
				box = (*iter)->box;
		}
	}
};

#endif /*CLOGICORCONTAINER_H_*/
