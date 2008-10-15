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

#include "geo_scene.h"
#include "geo_utils.h"

void Scene_createmeta(lua_State *L)
{
  struct luaL_reg methods[] = {
    {"add", Scene_add_object },
    {"__gc",Scene_destroy },
    {NULL, NULL},
  };
  luageo_createmeta(L, LUAGEO_SCENE, methods);
  lua_pop (L, 1);
}

int Scene_create(lua_State *L)
{
  Scene** sceneptr = (Scene **)lua_newuserdata(L, sizeof(Scene*));
  (*sceneptr) = new Scene();
  luageo_setmeta(L, LUAGEO_SCENE);
  return 1;
}

int Scene_destroy(lua_State *L)
{
  Scene** sceneptr = (Scene **)luaL_checkudata(L, 1, LUAGEO_SCENE);
  for (vector<CObject*>::iterator iter = (*sceneptr)->objects.begin(); 
       iter != (*sceneptr)->objects.end(); iter++)
    delete *iter;
  delete (*sceneptr);
}

int Scene_add_object(lua_State *L)
{
  Scene** sceneptr = (Scene **)luaL_checkudata(L, 1, LUAGEO_SCENE);
  luaL_checktype(L, 2, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  geo_obj* objptr = (geo_obj*)luaL_checkudata(L, -1, LUAGEO_OBJECT);
  lua_pop(L,1);
  CObject* cobj = objptr->object->clone();
  (*sceneptr)->objects.push_back(cobj);
  double depth = 0;
  geo_getdouble(L, "depth", &depth );
  cobj->dDepth = depth;
  return 1;
}
