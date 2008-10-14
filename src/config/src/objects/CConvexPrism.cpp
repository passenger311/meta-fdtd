
#include "CConvexPrism.h"

CConvexPrism::CConvexPrism() {} 

CConvexPrism::CConvexPrism(vector<vec2>& pointlist, double height)
{
  m_vPointlist.clear();
  for (vector<vec2>::iterator pt = pointlist.begin(); pt != pointlist.end(); pt++)
    m_vPointlist.push_back(*pt);
  m_dHalfHeight = height / 2;
}

CObject* CConvexPrism::clone() {
  return new CConvexPrism(m_vPointlist, m_dHalfHeight*2.0);
}

void CConvexPrism::preProcess()
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
	

bool CConvexPrism::PointInside(vec3& point)
{
  for (vector<line2>::iterator iline = m_lBorders.begin(); iline != m_lBorders.end(); iline++)
    {
      vec2 v2(point,VZ);
      if (iline->signedDistance(v2) < 0)
	return false;
    }
  return true;
}
