#include "CLogicOrContainer.h"

CLogicOrContainer::CLogicOrContainer() {

}

CLogicOrContainer::CLogicOrContainer(vector<CObject*>& objectlist)
{
  for (vector<CObject*>::iterator iter = objectlist.begin(); iter != objectlist.end(); iter++)
    m_vObjectList.push_back(*iter);
}

CObject* CLogicOrContainer::clone() {
  CLogicOrContainer* obj = new CLogicOrContainer();
  for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++)
    obj->m_vObjectList.push_back(*iter);
  return obj;
}

int CLogicOrContainer::Inside(pointframe& pfrm)
{
  if (!box.frameIntersection(pfrm))
    return 0;
  int insd = 0;
  for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++)
    insd = MAX(insd, (*iter)->Inside(pfrm));
  return insd;
}
	

bool CLogicOrContainer::PointInside(vec3& point)
{
  for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++)
    if ((*iter)->PointInside(point)) return true;
  return false;
}

bool CLogicOrContainer::addSubObject(CObject* object)
{ 
  m_vObjectList.push_back(object);
  return true;
}

void CLogicOrContainer::preProcess()
{
  for (vector<CObject*>::iterator iter = m_vObjectList.begin(); iter != m_vObjectList.end(); iter++) {
    (*iter)->preProcess();
    if (iter != m_vObjectList.begin())
      box = box.combineWith((*iter)->box);
    else
      box = (*iter)->box;
  }
}
