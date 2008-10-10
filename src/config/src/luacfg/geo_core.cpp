
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
#include "geo_grid.h"
#include "geo_objects.h"

extern "C" LUAGEO_API int luaopen_geo_core(lua_State *L);

LUAGEO_API int luaopen_geo_core(lua_State *L)
{
  struct luaL_reg reg[] = {
    {"Scene_create", Scene_create},
    {"Scene_destroy", Scene_destroy},
    {"Grid_create", Grid_create},
    {"Grid_destroy", Grid_destroy},
    {"Sphere_create", CSphere_create},
    {"Object_destroy", CObject_destroy},
    {NULL, NULL},
  };

  Scene_create_metatable(L);
  Grid_create_metatable(L);
  Objects_create_metatable(L);

  luaL_openlib (L, "geo_core", reg, 0);
  luageo_set_info (L);
  return 1;
}


