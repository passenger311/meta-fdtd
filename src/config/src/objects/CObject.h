#ifndef COBJECT_H_
#define COBJECT_H_

#define DEBUG

#include "algebra3extension.h"
#include <string>

using namespace std;

class CObject
{
public:
	CObject();
	virtual ~CObject();
	virtual CObject* clone() { return new CObject(); }
	
public:
	frame box;
	int iLastInside;
	string name;
	// [AH] --->
	double fWeight;
	double dDepth;
	// <---

	virtual void preProcess();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
	virtual bool InsideFrame(vec3& point);
	virtual bool InsideFrame(pointframe& pfrm);
	virtual bool addSubObject(CObject* object) { return false; }
};



struct PredicateCObjectDepth
{
     bool operator()(CObject* const& one, CObject* const& two )
     {
          return one->dDepth < two->dDepth;
     }
};


#endif /*COBJECT_H_*/
