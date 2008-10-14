#include "CSimpleTransform.h"

CSimpleTransform::CSimpleTransform()
{
  oObject = NULL;
  vTranslation = vec3(0,0,0);
  vRotationAxis = vec3(1,0,0);
  dRotationAngle = 0;
  dScaleFactor = 1;
}

CObject* CSimpleTransform::clone() {
  return new CSimpleTransform(oObject->clone(),vTranslation,vRotationAxis,dRotationAngle,dScaleFactor); 
}

CSimpleTransform::CSimpleTransform(CObject* object, vec3& translation, vec3& rotationaxis, double angle, double scale)
{
  oObject = object;
  vTranslation = translation;
  vRotationAxis = rotationaxis;
  dRotationAngle = angle;
  dScaleFactor = scale;
}

bool CSimpleTransform::addSubObject(CObject* object)
{ 
  if (oObject == NULL)
    oObject = object;
  else
    return false;
  return true;
}

void CSimpleTransform::preProcess()
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

bool CSimpleTransform::PointInside(vec3& point)
{
  vec3 pt = m_mTransform * point;
  return oObject->PointInside(pt);
}


int CSimpleTransform::Inside(pointframe& pfrm)
{
	if (!box.frameIntersection(pfrm))
		return 0;
	pointframe pff;
	pff = pfrm;
	pff.transformByMatrix(m_mTransform,dScaleFactor); //TODO: 1/dScaleFactor??
	return oObject->Inside(pff);
}
