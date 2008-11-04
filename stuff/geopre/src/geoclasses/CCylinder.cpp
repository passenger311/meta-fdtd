#include "CCylinder.h"

inline bool CCylinder::PointInside(vec3& point)
{
	return point[VX]*point[VX]+point[VY]*point[VY] <= m_dRadiusSquare &&
		-m_dHalfHeight <= point[VZ] && m_dHalfHeight >= point[VZ];
}
