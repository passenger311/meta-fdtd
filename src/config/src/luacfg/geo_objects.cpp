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
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CSphere*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CBezierPrism_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CBezierPrism*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}


int CBox_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CBox*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CConvexPrism_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CConvexPrism*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CCylinder_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CCylinder*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CEllipsoid_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CEllipsoid*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicAndNotObject_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CLogicAndNotObject*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicAndObject_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CLogicAndObject*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicOrContainer_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CLogicOrContainer*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicOrObject_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CLogicOrObject*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicXOrObject_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CLogicXOrObject*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimplePrism_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CSimplePrism*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimpleRotationZ_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CSimpleRotationZ*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimpleTransform_create(lua_State *L) {
  CObject** objptr = (CObject **)lua_newuserdata(L, sizeof(CSimpleTransform*));
  //  (*objptr) = new CSphere();
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}
