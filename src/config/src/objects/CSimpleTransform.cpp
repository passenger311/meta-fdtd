#include "CSimpleTransform.h"

CSimpleTransform::~CSimpleTransform()
{
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
