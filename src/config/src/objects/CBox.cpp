#include "CBox.h"

CBox::CBox(vec3& pos, vec3& size) {
  vec3 sv = pos - 0.5 * size;
  vec3 ev = pos + 0.5 * size;
  box = frame(sv,ev);
#ifdef DEBUG
  printf("Box: from = %f %f %f to = %f %f %f\n",sv[0],sv[1],sv[2],ev[0],ev[1],ev[2]);
#endif
}

CBox::CBox(frame& f) {
  box = f;
}

CObject* CBox::clone() { return new CBox(box); }
	
bool CBox::PointInside(vec3& point)
{
  return box.pointInFrame(point);
}
