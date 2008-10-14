#ifndef CLOGICORCONTAINER_H_
#define CLOGICORCONTAINER_H_

#include "CObject.h"
#include <vector>

using namespace std;

class CLogicOrContainer : public CObject
{
protected:
	vector<CObject*> m_vObjectList;	

public:
	CLogicOrContainer();
	CLogicOrContainer(vector<CObject*>& objectlist);
	virtual CObject* clone();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
	virtual bool addSubObject(CObject* object);
	virtual void preProcess();
};

#endif /*CLOGICORCONTAINER_H_*/
