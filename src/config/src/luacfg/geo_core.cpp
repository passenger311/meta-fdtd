
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

extern "C" LUAGEO_API int luaopen_geo(lua_State *L);

LUAGEO_API int luaopen_geo(lua_State *L)
{
  struct luaL_reg reg[] = {
    {"Scene", Scene_create},
    {"Grid", Grid_create},
    // {"FileVTK", FileVTK_create},
    // {"FileIN", FileIN_create},
    {"BezierPrism", CBezierPrism_create},
    {"Box", CBox_create },
    {"ConvexPrism", CConvexPrism_create },
    {"Cylinder", CCylinder_create },
    {"Ellipsoid", CEllipsoid_create },
    {"BinaryAndNot", CLogicAndNotObject_create },
    {"BinaryAnd", CLogicAndObject_create },
    {"CollectionOr", CLogicOrContainer_create },
    {"BinaryOr", CLogicOrObject_create },
    {"BinaryXOr", CLogicXOrObject_create },
    {"Prism", CSimplePrism_create },
    {"RotationZ", CSimpleRotationZ_create },
    {"Transform", CSimpleTransform_create },
    {"Sphere", CSphere_create },
    {NULL, NULL},
  };

  Scene_createmeta(L);
  Grid_createmeta(L);
  // File_createmeta(L);
  Objects_createmeta(L);

  luaL_openlib (L, "geo", reg, 0);
  luageo_set_info (L);
  return 1;
}


