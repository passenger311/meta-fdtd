#include "CBox.h"

CBox::CBox(vec3& pos, vec3& size) {
  vec3 sv = pos - 0.5 * size;
  vec3 ev = pos + 0.5 * size;
  box = frame(sv,ev);
}

CBox::CBox(frame& f) {
  box = f;
}

CObject* CBox::clone() { return new CBox(box); }
	
bool CBox::PointInside(vec3& point)
{
  return box.pointInFrame(point);
}
