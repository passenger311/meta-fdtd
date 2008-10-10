#ifndef _GEO_SCENE_
#define _GEO_SCENE_
/* ----------------------------------------------------------------------- */

#include "../Scene.h"
#include "../objects/CObject.h"

void Scene_create_metatable(lua_State *L);
int Scene_create(lua_State *L);
int Scene_destroy(lua_State *L);
int Scene_add_object(lua_State *L);

/* ----------------------------------------------------------------------- */
#endif
