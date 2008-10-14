#include "CBinaryObject.h"

CBinaryObject::CBinaryObject()
{
	oObject1 = NULL;
	oObject2 = NULL;
}

CBinaryObject::CBinaryObject(CObject* obj1, CObject* obj2)
{
	oObject1 = obj1;	
	oObject2 = obj2;	
}
