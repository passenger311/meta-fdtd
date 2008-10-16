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

#include "geo_utils.h"
#include "geo_file.h"


void File_createmeta(lua_State *L)
{
    struct luaL_reg methods[] = {
        { "__gc", File_destroy },
	{NULL, NULL},
    };
    luageo_createmeta(L, LUAGEO_FILE, methods);
    lua_pop (L, 1);
}

int FileVTK_create(lua_State *L)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  const char* name = luaL_checkstring(L, -1);
  lua_pop(L,1); 
  FileHandler** fileptr = (FileHandler**)lua_newuserdata(L, sizeof(void*));
  (*fileptr) = new FileHandlerVTK(name);
  luageo_setmeta(L, LUAGEO_FILE);
  return 1;
}

int FileIN_create(lua_State *L)
{
  luaL_checktype(L, 1, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  const char* name = luaL_checkstring(L, -1);
  lua_pop(L,1);
  double comps;
  geo_getdouble(L, "comps", &comps );
  FileHandler** fileptr = (FileHandler**)lua_newuserdata(L, sizeof(void*));
  (*fileptr) = new FileHandlerFortranIN(name,(int)(comps+0.5));
  luageo_setmeta(L, LUAGEO_FILE);
  return 1;
}

int File_destroy(lua_State *L)
{
  FileHandler** fileptr = (FileHandler **)luaL_checkudata(L, 1, LUAGEO_FILE);
  delete (*fileptr);
}
