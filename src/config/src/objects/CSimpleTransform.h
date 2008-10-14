#ifndef CSIMPLETRANSFORM_H_
#define CSIMPLETRANSFORM_H_

#include "CObject.h"

class CSimpleTransform : public CObject
{
protected:
	mat4 m_mTransform;
	mat4 m_mTransformInverse;
public:
	vec3 vTranslation;
	vec3 vRotationAxis;
	double dRotationAngle;
	double dScaleFactor;
	CObject* oObject;
	
public:
	CSimpleTransform();
	CSimpleTransform(CObject* object, vec3& translation, vec3& rotationaxis, double angle, double scale);
	virtual CObject* clone(); 
	virtual bool addSubObject(CObject* object);
	virtual void preProcess();
	virtual int Inside(pointframe& pfrm);
	virtual bool PointInside(vec3& point);
};

#endif /*CSIMPLETRANSFORM_H_*/
