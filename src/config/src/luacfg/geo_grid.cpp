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
#include "geo_grid.h"



int Grid_destroy(lua_State *L)
{
  Grid** gridptr = (Grid **)luaL_checkudata(L, 1, LUAGEO_GRID);
  delete (*gridptr);
}

int Grid_write_output(lua_State *L)
{
  Grid** gridptr = (Grid **)luaL_checkudata(L, 1, LUAGEO_GRID);
  Scene** sceneptr = (Scene **)luaL_checkudata(L, 2, LUAGEO_SCENE);
  FileHandler** fileptr = (FileHandler**)luaL_checkudata(L, 2, LUAGEO_FILE);
  (*gridptr)->generateOutput(*sceneptr,*fileptr);
  return 1;
}

void Grid_createmeta(lua_State *L)
{
    struct luaL_reg methods[] = {
        {"write", Grid_write_output },
	{"__gc", Grid_destroy },
	{NULL, NULL},
    };
    luageo_createmeta(L, LUAGEO_GRID, methods);
    lua_pop (L, 1);
}

int Grid_create(lua_State *L)
{
  double pos0[3] = { 0.,0.,0. };
  double pos1[3] = { 1.,1.,1. };
  double cells[3] = { 10.,10.,10. };
  if ( lua_istable(L,1) ) {
    geo_getvec3(L, "from", pos0);
    geo_getvec3(L, "to", pos1);
    geo_getvec3(L, "cells", cells );
  }  
  vec3 from(pos0[0],pos0[1],pos0[2]);
  vec3 to(pos1[0],pos1[1],pos1[2]);
  frame f(from,to);
  Grid** gridptr = (Grid **)lua_newuserdata(L, sizeof(Grid*));
  (*gridptr) = new Grid(f,(int)cells[0],(int)cells[1],(int)cells[2]);
  luageo_setmeta(L, LUAGEO_GRID);
  return 1;
}
