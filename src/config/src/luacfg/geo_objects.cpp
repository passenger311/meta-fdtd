extern "C" {

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "lua.h"
#include "lauxlib.h"
#if ! defined (LUA_VERSION_NUM) || LUA_VERSION_NUM < 501
#include "compat-5.1.h"
#endif

#include "luageo.h"

}

#include "geo_objects.h"

void Objects_create_metatable(lua_State *L)
{
    struct luaL_reg methods[] = {
	{NULL, NULL},
    };

    luageo_createmeta(L, LUAGEO_SCENE, methods);
    lua_pop (L, 1);
}


int CObject_destroy(lua_State *L)
{
  CObject** objptr = (CObject **)luaL_checkudata(L, 1, LUAGEO_OBJECT);
  delete (*objptr);
}


int CSphere_create(lua_State *L) {

}

int CBezierPrism_create(lua_State *L) {

}

int CBinaryObject_create(lua_State *L) {

}

int CBox_create(lua_State *L) {

}

int CConvexPrism_create(lua_State *L) {

}

int CCylinder_create(lua_State *L) {

}

int CEllipsoid_create(lua_State *L) {

}

int CLogicAndNotObject_create(lua_State *L) {

}

int CLogicAndObject_create(lua_State *L) {

}

int CLogicOrContainer_create(lua_State *L) {

}

int CLogicOrObject_create(lua_State *L) {

}

int CLogicXOrObject_create(lua_State *L) {

}

int CSimplePrism_create(lua_State *L) {

}

int CSimpleRotationZ_create(lua_State *L) {

}

int CSimpleTransform_create(lua_State *L) {

}
