
#ifndef _LUAGEO_UTILS_
#define _LUAGEO_UTILS_

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
#include "../objects/CObject.h"

using namespace std;

struct geo_obj_s;

typedef struct geo_obj_s {
  short  refcount;
  CObject* object;
  short numchild;
  CObject** childs;
} geo_obj;

int geo_getvec3(lua_State *L, const char* name, double val[3]);
int geo_getdouble(lua_State *L, const char* name, double* val);
int geo_getbool(lua_State *L, const char* name, bool* val);
int geo_getunary(lua_State *L, CObject** obj);
int geo_getbinary(lua_State *L, CObject** obj1, CObject** obj2);
int geo_getcollection(lua_State *L, vector<CObject*>& vec);
int geo_getpointlist2(lua_State *L, const char* name, vector<vec2>& points2);

#endif
