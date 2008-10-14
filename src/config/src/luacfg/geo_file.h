#ifndef _GEO_FILE_
#define _GEO_FILE_
/* ----------------------------------------------------------------------- */

#include "../File.h"
#include "../Scene.h"

void File_createmeta(lua_State *L);

int FileVTK_create(lua_State *L);
int FileIN_create(lua_State *L);
int File_destroy(lua_State *L);

/* ----------------------------------------------------------------------- */
#endif
