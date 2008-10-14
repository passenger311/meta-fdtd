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

#include <vector>
#include "geo_objects.h"

void Objects_createmeta(lua_State *L)
{
    struct luaL_reg methods[] = {
      { "__gc", geo_obj_collect },
      {NULL, NULL},
    };

    luageo_createmeta(L, LUAGEO_OBJECT, methods);
    lua_pop (L, 1);
}


void geo_obj_addchild(geo_obj* obj, CObject* child) {
  obj->numchild++;
  CObject** new_childs = new CObject*[obj->numchild];
  for(int i=0; i<obj->numchild-1; i++) { new_childs[i] = obj->childs[i]; }
  new_childs[obj->numchild-1] = child;
  delete[] obj->childs;
  obj->childs = new_childs;
}

int geo_obj_collect(lua_State *L)
{
  geo_obj* obj = (geo_obj*)luaL_checkudata(L, 1, LUAGEO_OBJECT);
  for(int i = 0; i<obj->numchild; i++ ) { delete obj->childs[i]; }
  if ( obj->childs ) delete[] obj->childs;
  delete obj->object; 
}

geo_obj* geo_obj_create(lua_State *L,CObject* cobject)
{
  geo_obj* obj = (geo_obj*)lua_newuserdata(L, sizeof(geo_obj));
  obj->refcount = 1;
  obj->object = cobject;
  obj->childs = NULL;
  obj->numchild = 0;
  return obj;
}

int CSphere_create(lua_State *L) {
  double pos[3] = { 0.,0.,0. };
  double rad = 1.;
  if ( lua_istable(L,1) ) {
    geo_getvec3(L, "at", pos);
    geo_getdouble(L, "radius", &rad );
  }
  vec3 p(pos[0],pos[1],pos[2]);  
  geo_obj_create(L, new CSphere(p,rad));
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CCylinder_create(lua_State *L) {
  double pos[3] = { 0.,0.,0. };
  double rad = 1.;
  double height = 1.;
  if ( lua_istable(L,1) ) {
    geo_getvec3(L, "at", pos);
    geo_getdouble(L, "height", &height);
    geo_getdouble(L, "radius", &rad );
  }  
  vec3 p(pos[0],pos[1],pos[2]);  
  geo_obj_create(L, new CCylinder(p,rad,height)); 
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CBox_create(lua_State *L) {
  double pos[3] = {0.,0.,0.};
  double size[3] = {1.,1.,1.};
  if ( lua_istable(L,1) ) {
    geo_getvec3(L, "at", pos);
    geo_getvec3(L, "size", size );
  }  
  vec3 p(pos[0],pos[1],pos[2]);
  vec3 s(size[0],size[1],size[2]);
  geo_obj_create(L, new CBox(p,s)); 
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CEllipsoid_create(lua_State *L) {
  double pos[3] = {0.,0.,0.};
  double size[3] = {1.,1.,1.};
  if ( lua_istable(L,1) ) {
    geo_getvec3(L, "at", pos);
    geo_getvec3(L, "size", size );
  }  
  vec3 p(pos[0],pos[1],pos[2]);
  geo_obj_create(L, new CEllipsoid(p,size[0],size[1],size[2])); 
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicOrObject_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  CObject* obj2;
  geo_getbinary(L,&obj1,&obj2);
  geo_obj* obj = geo_obj_create(L, new CLogicOrObject(obj1,obj2));
  geo_obj_addchild(obj,obj1);
  geo_obj_addchild(obj,obj2);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicAndObject_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  CObject* obj2;
  geo_getbinary(L,&obj1,&obj2);
  geo_obj* obj = geo_obj_create(L, new CLogicAndObject(obj1,obj2));
  geo_obj_addchild(obj,obj1);
  geo_obj_addchild(obj,obj2);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicAndNotObject_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  CObject* obj2;
  geo_getbinary(L,&obj1,&obj2);
  geo_obj* obj = geo_obj_create(L, new CLogicAndNotObject(obj1,obj2));
  geo_obj_addchild(obj,obj1);
  geo_obj_addchild(obj,obj2);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}


int CLogicXOrObject_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  CObject* obj2;
  geo_getbinary(L,&obj1,&obj2);
  geo_obj* obj = geo_obj_create(L, new CLogicXOrObject(obj1,obj2));
  geo_obj_addchild(obj,obj1);
  geo_obj_addchild(obj,obj2);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CLogicOrContainer_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  vector<CObject*> objv;
  geo_getcollection(L,objv);
  geo_obj* obj = geo_obj_create(L, new CLogicOrContainer(objv));
  for (vector<CObject*>::iterator iter = objv.begin(); iter != objv.end(); iter++)
    geo_obj_addchild(obj,*iter);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimpleRotationZ_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  geo_getunary(L,&obj1);
  geo_obj* obj = geo_obj_create(L, new CSimpleRotationZ(obj1));
  geo_obj_addchild(obj,obj1);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimpleTransform_create(lua_State *L) {
  double move[3] = { 0.,0.,0. };
  double axis[3] = { 0.,0.,1. };
  double angle = 0.; 
  double scale = 1.;
  luaL_checktype(L, 1, LUA_TTABLE);
  CObject* obj1;
  geo_getunary(L,&obj1);
  geo_getvec3(L, "move", move);
  geo_getvec3(L, "axis", axis);
  geo_getdouble(L, "angle", &angle);
  geo_getdouble(L, "scale", &scale);
  vec3 vmove(move[0],move[1],move[2]);
  vec3 vaxis(axis[0],axis[1],axis[2]);
  geo_obj* obj = geo_obj_create(L, new CSimpleTransform(obj1,vmove,vaxis,angle,scale));
  geo_obj_addchild(obj,obj1);
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CConvexPrism_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  vector<vec2> points2;
  double height = 1.;
  geo_getpointlist2(L, "points", points2);
  geo_getdouble(L, "height", &height);
  geo_obj_create(L, new CConvexPrism(points2,height));
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CSimplePrism_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  vector<vec2> points2;
  double height = 1.;
  geo_getpointlist2(L, "points", points2);
  geo_getdouble(L, "height", &height);
  geo_obj_create(L, new CSimplePrism(points2,height));
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}

int CBezierPrism_create(lua_State *L) {
  luaL_checktype(L, 1, LUA_TTABLE);
  vector<vec2> points2;
  double height = 1.;
  double steps = 10.;
  geo_getpointlist2(L, "points", points2);
  geo_getdouble(L, "height", &height);
  geo_getdouble(L, "steps", &steps);
  geo_obj_create(L, new CBezierPrism(points2,height,(int)steps));
  luageo_setmeta(L, LUAGEO_OBJECT);
  return 1;
}



