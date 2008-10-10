#ifndef _GEO_OBJECTS_
#define _GEO_OBJECTS_
/* ----------------------------------------------------------------------- */

#include "../objects/allobjects.h"

void Objects_create_metatable(lua_State *L);

int CObject_destroy(lua_State *L);

int CSphere_create(lua_State *L);
int CBezierPrism_create(lua_State *L);
int CBinaryObject_create(lua_State *L);
int CBox_create(lua_State *L);
int CConvexPrism_create(lua_State *L);
int CCylinder_create(lua_State *L);
int CEllipsoid_create(lua_State *L);
int CLogicAndNotObject_create(lua_State *L);
int CLogicAndObject_create(lua_State *L);
int CLogicOrContainer_create(lua_State *L);
int CLogicOrObject_create(lua_State *L);
int CLogicXOrObject_create(lua_State *L);
int CObject_create(lua_State *L);
int CSimplePrism_create(lua_State *L);
int CSimpleRotationZ_create(lua_State *L);
int CSimpleTransform_create(lua_State *L);
int CSphere_create(lua_State *L);

/* ----------------------------------------------------------------------- */
#endif
