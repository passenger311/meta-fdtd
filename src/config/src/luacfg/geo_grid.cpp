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
  luaL_checktype(L, 2, LUA_TTABLE);
  lua_pushnumber(L,1);
  lua_gettable(L,-2);
  FileHandler** fileptr = (FileHandler**)luaL_checkudata(L, -1, LUAGEO_FILE);
  lua_pop(L,1);
  lua_pushnumber(L,2);
  lua_gettable(L,-2);
  Scene** sceneptr = (Scene **)luaL_checkudata(L, -1, LUAGEO_SCENE);
  lua_pop(L,1);
  char* method = "default";
  geo_getstring(L, "method", &method );
  bool silent = true;
  geo_getbool(L, "silent", &silent );
  (*gridptr)->generateOutput(*sceneptr,*fileptr,method,silent);
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
  double from[3] = { -1.,-1.,-1. };
  double to[3] = { 1.,1.,1. };
  double cells[3] = { -1,-1,-1 };
  double cagescale = 2.;
  double cagediv[3] = { 1.,1.,1. };
  double subdiv[3] = { 5.,5.,5. };
  bool cagefps = true;
  bool nosub = false;
  bool forcesub = false;
  bool yee = true;
  luaL_checktype(L, 1, LUA_TTABLE);
  geo_getvec3(L, "from", from);
  geo_getvec3(L, "to", to);
  geo_getvec3(L, "cells", cells );
  geo_getdouble(L, "cagescale", &cagescale );
  geo_getvec3(L, "cagediv", cagediv );
  geo_getbool(L, "cagefp", &cagefps );
  geo_getvec3(L, "subdiv", subdiv );
  geo_getbool(L, "nosub", &nosub );
  geo_getbool(L, "forcesub", &forcesub );
  geo_getbool(L, "yee", &yee );   
  vec3 vfrom(from[0],from[1],from[2]);
  vec3 vto(to[0],to[1],to[2]);
  frame fr(vfrom,vto);
  Grid** gridptr = (Grid **)lua_newuserdata(L, sizeof(Grid*));
  (*gridptr) = new Grid(fr,(int)cells[0],(int)cells[1],(int)cells[2]);
  (*gridptr)->iSubGriddingDivX = (int)subdiv[0];
  (*gridptr)->iSubGriddingDivY = (int)subdiv[1];
  (*gridptr)->iSubGriddingDivZ = (int)subdiv[2];
  (*gridptr)->iPointframeDivisionsX = (int)cagediv[0];
  (*gridptr)->iPointframeDivisionsY = (int)cagediv[1];
  (*gridptr)->iPointframeDivisionsZ = (int)cagediv[2];
  (*gridptr)->dPointframeScale = cagescale;
  (*gridptr)->bPointframeFacepoints = cagefps;
  (*gridptr)->bAlwaysSubgridding = forcesub;
  (*gridptr)->bNoSubgridding = nosub;
  (*gridptr)->bYeeGrid = yee;
  luageo_setmeta(L, LUAGEO_GRID);
  return 1;
}
