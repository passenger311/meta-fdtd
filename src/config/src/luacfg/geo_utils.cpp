
#include "geo_utils.h"

/* --- Get a double vec3 from a table */
int geo_getvec3(lua_State *L, const char* name, double val[3]) {
  int ok = 0;
  int i;
  lua_getfield(L, -1, name);
  if ( lua_istable(L,-1) ) { 
    ok = 1;
    for(i=0; i<3; i++) {
      lua_pushnumber(L,i+1);
      lua_gettable(L,-2);
      if ( lua_isnumber(L,-1) ) { val[i] = lua_tonumber(L,-1); } 
      lua_pop(L,1); /* pop number */
    } 
  }
  lua_pop(L,1); /* pop field */
  return ok;
}


/* --- Get a double from a table */
int geo_getdouble(lua_State *L, const char* name, double* val) {
  int ok = 0;
  luaL_checktype(L, -1, LUA_TTABLE);
  lua_getfield(L, -1, name);
  if ( lua_isnumber(L,-1) ) { 
    ok = 1;
    *val = lua_tonumber(L,-1);
  }
  lua_pop(L,1); /* pop field */
  return ok;
}

/* --- Get a string from a table */
int geo_getstring(lua_State *L, const char* name, char** val) {
  int ok = 0;
  luaL_checktype(L, -1, LUA_TTABLE);
  lua_getfield(L, -1, name);
  if ( lua_isstring(L,-1) ) { 
    ok = 1;
    *val = (char*) lua_tostring(L,-1);
  }
  lua_pop(L,1); /* pop field */
  return ok;
}


/* --- Get a bool from a table */
int geo_getbool(lua_State *L, const char* name, bool* val) {
  int ok = 0;
  luaL_checktype(L, -1, LUA_TTABLE);
  lua_getfield(L, -1, name);
  if ( lua_isboolean(L,-1) ) { 
    ok = 1;
    *val = lua_toboolean(L,-1) == 1 ? true : false;
  }
  lua_pop(L,1); /* pop field */
  return ok;
}



/* --- Get single object from table */
int geo_getunary(lua_State *L, CObject** obj) {
  luaL_checktype(L, -1, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  *obj = ((geo_obj*)luaL_checkudata(L, -1, LUAGEO_OBJECT))->object->clone();
  lua_pop(L,1);
}


/* --- Get pair of objects from table */
int geo_getbinary(lua_State *L, CObject** obj1, CObject** obj2) {
  luaL_checktype(L, -1, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  *obj1 = ((geo_obj*)luaL_checkudata(L, -1, LUAGEO_OBJECT))->object->clone();
  lua_pop(L,1);
  lua_pushnumber(L,2);
  lua_gettable(L,-2);
  *obj2 = ((geo_obj*)luaL_checkudata(L, -1, LUAGEO_OBJECT))->object->clone();
  lua_pop(L,1);
}

/* get table size */
int geo_getn(lua_State *L, int idx) {
  lua_pushstring(L, "table"); 
  lua_gettable(L, LUA_GLOBALSINDEX);
  lua_pushstring(L, "getn");
  lua_gettable(L, -2);
  lua_pushvalue(L, idx);  
  lua_call(L, 1, 1);
  int n = (int)lua_tonumber(L,-1);
  lua_pop(L,2);
  return n;
}

/* --- Get a list of objects from table */
int geo_getcollection(lua_State *L, vector<CObject*>& vec) {
  luaL_checktype(L, -1, LUA_TTABLE);
  int n = geo_getn(L,1);
  CObject* obj;
  for (int i = 1; i <= n; i++) {
    lua_pushnumber(L,i);
    lua_gettable(L,-2);
    obj = ((geo_obj*)luaL_checkudata(L, -1, LUAGEO_OBJECT))->object->clone();
    vec.push_back(obj);
    lua_pop(L,1);
  }
}

/* --- Read a list of 2d points from field */
int  geo_getpointlist2(lua_State *L, const char* name, vector<vec2>& points2) {
  luaL_checktype(L, -1, LUA_TTABLE);
  int ok = 0;
  int i;
  lua_getfield(L, 1, name);
  if ( lua_istable(L,-1) ) { 
    int n = geo_getn(L,2);
    for (int i = 1; i <= n; i++) {
      lua_pushnumber(L,i);
      lua_gettable(L,-2);
      lua_pushnumber(L,1);
      lua_gettable(L,-2);
      double x = lua_tonumber(L,-1);
      lua_pop(L,1);
      lua_pushnumber(L,2);
      lua_gettable(L,-2);
      double y = lua_tonumber(L,-1);
      lua_pop(L,1);
      vec2 point(x,y);
      points2.push_back(point);
      lua_pop(L,1);
    }
    ok = 1;
  }
  lua_pop(L,1); /* pop field */
  return ok;
}
