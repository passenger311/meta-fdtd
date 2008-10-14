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

#include "geo_grid.h"


int Grid_write_output(lua_State *L)
{
	Grid** gridptr = (Grid **)luaL_checkudata(L, 1, LUAGEO_GRID);
	//	CObject** objectptr = (CObject **)luaL_checkudata(L, 2, LUAGEO_OBJECT);
	//	luaL_argcheck(L, grid != NULL, 1, LUAGEO_PREFIX"grid expected");
	

	return 1;
}


void Grid_createmeta(lua_State *L)
{
    struct luaL_reg methods[] = {
        {"write", Grid_write_output },
	{NULL, NULL},
    };

    luageo_createmeta(L, LUAGEO_GRID, methods);
    lua_pop (L, 1);
}

int Grid_create(lua_State *L)
{
	Grid** gridptr = (Grid **)lua_newuserdata(L, sizeof(Grid*));
	(*gridptr) = new Grid();
	luageo_setmeta(L, LUAGEO_GRID);
	return 1;
}

int Grid_destroy(lua_State *L)
{
  Grid** gridptr = (Grid **)luaL_checkudata(L, 1, LUAGEO_GRID);
  delete (*gridptr);
}
