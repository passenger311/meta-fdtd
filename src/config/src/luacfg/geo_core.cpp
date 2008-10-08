/*

*/

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

#include "../Scene.h"

#define LUAGEO_SCENE "Scene object"
#define LUAGEO_GRID "Grid object"

 
extern "C" LUAGEO_API int luaopen_geo(lua_State *L);

/*
** close environment object.
*/
static int scene_hello(lua_State *L)
{
	Scene *scene = (Scene *)luaL_checkudata(L, 1, LUAGEO_SCENE);
	luaL_argcheck(L, scene != NULL, 1, LUAGEO_PREFIX"scene expected");
	lua_pushboolean(L, 1);
	return 1;
}

/*
** create metatable for scene
*/
static void create_metatable_scene(lua_State *L)
{
    struct luaL_reg methods[] = {
        {"hello", scene_hello},
	{NULL, NULL},
    };

    luageo_createmeta(L, LUAGEO_SCENE, methods);
    lua_pop (L, 1);
}

/*
** create scene and return it
*/
static int create_new_scene (lua_State *L)
{
	Scene *scene = (Scene *)lua_newuserdata(L, sizeof(Scene));
	luageo_setmeta(L, LUAGEO_SCENE);
	return 1;
}

/*
** register all objects 
*/

  LUAGEO_API int luaopen_geo(lua_State *L)
  {
    struct luaL_reg reg[] = {
      {"Scene", create_new_scene},
      {NULL, NULL},
    };
    create_metatable_scene(L);
    luaL_openlib (L, "geo", reg, 0);
    luageo_set_info (L);
    return 1;
  }


