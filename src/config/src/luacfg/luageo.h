
#ifndef _LUAGEO_
#define _LUAGEO_

#ifndef LUAGEO_API
#define LUAGEO_API
#endif

#define LUAGEO_PREFIX "Geo:"

#define LUAGEO_SCENE "Geo:Scene"
#define LUAGEO_OBJECT "Geo:Object"
#define LUAGEO_GRID "Geo:Grid"

LUAGEO_API int luageo_error (lua_State *L, const char *err);
LUAGEO_API int luageo_createmeta (lua_State *L, const char *name, const luaL_reg *methods);
LUAGEO_API void luageo_setmeta (lua_State *L, const char *name);
LUAGEO_API void luageo_set_info (lua_State *L);

typedef struct {
        short  closed;
} pseudo_data;



#endif
