#ifndef _GEO_GRID_
#define _GEO_GRID_
/* ----------------------------------------------------------------------- */

#include "../Grid.h"
#include "../Scene.h"

void Grid_create_metatable(lua_State *L);

int Grid_create(lua_State *L);
int Grid_destroy(lua_State *L);

/* ----------------------------------------------------------------------- */
#endif
