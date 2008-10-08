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
	CSimpleTransform()
	{
		oObject = NULL;
		vTranslation = vec3(0,0,0);
		vRotationAxis = vec3(1,0,0);
		dRotationAngle = 0;
		dScaleFactor = 1;
	}
	CSimpleTransform(CObject* object, vec3& translation, vec3& rotationaxis, double angle, double scale)
	{
		oObject = object;
		vTranslation = translation;
		vRotationAxis = rotationaxis;
		dRotationAngle = angle;
		dScaleFactor = scale;
	}
	virtual ~CSimpleTransform();

	virtual bool addSubObject(CObject* object)
	{ 
		if (oObject == NULL)
			oObject = object;
		else
			return false;
		return true;
	}
	
	virtual void preProcess()
	{
		oObject->preProcess();
		m_mTransformInverse = translation3D(vTranslation) 
				* rotation3D(vRotationAxis,dRotationAngle) 
				* scaling3D(vec3(dScaleFactor,dScaleFactor,dScaleFactor));
		m_mTransform = m_mTransformInverse.inverse();
		vec3* pts = oObject->box.getEdgePoints();
		for (int i=0; i<8 ; i++)
			pts[i] = m_mTransformInverse * pts[i];
		box.getFromPointlist(pts,8);
		delete [] pts;
	}
	virtual int Inside(pointframe& pfrm);
	
	virtual bool PointInside(vec3& point)
	{
		vec3 pt = m_mTransform * point;
		return oObject->PointInside(pt);
	}
};

#endif /*CSIMPLETRANSFORM_H_*/
