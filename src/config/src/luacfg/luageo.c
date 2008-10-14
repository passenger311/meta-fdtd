
#include <string.h>
  
#include "lua.h"
#include "lauxlib.h"
#if ! defined (LUA_VERSION_NUM) || LUA_VERSION_NUM < 501
#include "compat-5.1.h"
#endif
  
#include "luageo.h"

/* --- Error */
LUAGEO_API int luageo_error(lua_State *L, const char *err) {
  lua_pushnil(L);
  lua_pushstring(L, err);
  return 2;
}


/* --- Return the name of the object's metatable; used by `tostring' */
static int luageo_tostring (lua_State *L) {
  char buff[100];
  pseudo_data *obj = (pseudo_data *)lua_touserdata (L, 1);
  if (obj->closed)
    strcpy (buff, "closed");
  else
    sprintf (buff, "%p", (void *)obj);
  lua_pushfstring (L, "%s (%s)", lua_tostring(L,lua_upvalueindex(1)), buff);
  return 1;
}


/* --- Create a metatable <name> and leave it on top of the stack */
LUAGEO_API int luageo_createmeta (lua_State *L, const char *name, const luaL_reg *methods) {
  if (!luaL_newmetatable (L, name))
    return 0;
  
  /* define methods */
  luaL_openlib (L, NULL, methods, 0);
  
  /* define metamethods */
  lua_pushliteral (L, "__gc");
  lua_pushcfunction (L, methods->func);
  lua_settable (L, -3);
  
  lua_pushliteral (L, "__index");
  lua_pushvalue (L, -2);
  lua_settable (L, -3);
  
  lua_pushliteral (L, "__tostring");
  lua_pushstring (L, name);
  lua_pushcclosure (L, luageo_tostring, 1);
  lua_settable (L, -3);
  
  lua_pushliteral (L, "__metatable");
  lua_pushliteral (L, LUAGEO_PREFIX"you're not allowed to get this metatable");
  lua_settable (L, -3);
  
  return 1;
}


/* --- Assign metatable <name> to the object on top of the stack */
LUAGEO_API void luageo_setmeta (lua_State *L, const char *name) {
  luaL_getmetatable (L, name);
  lua_setmetatable (L, -2);
}

/* --- Assumes the table is on top of the stack */
LUAGEO_API void luageo_set_info (lua_State *L) {
  lua_pushliteral (L, "_COPYRIGHT");
  lua_pushliteral (L, "J Hamm");
  lua_settable (L, -3);
  lua_pushliteral (L, "_DESCRIPTION");
  lua_pushliteral (L, "GEO is a preprocessor for FDTD codes");
  lua_settable (L, -3);
  lua_pushliteral (L, "_VERSION");
  lua_pushliteral (L, "GEO 1.0.0");
  lua_settable (L, -3);
}


