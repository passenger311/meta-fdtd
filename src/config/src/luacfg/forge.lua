local INF = {
------------------------------------------------------------------------------
------------------------------------------------------------------------------
file        = "forge",
version     = "0.9",
date        = "9/1/2008",
author      = "j.hamm",
license     = "x11", 
lua_ver     = "5.1"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
}

------------------------------------------------------------------------------

local PREFIX   = "LUAOPEN_API"
local PATH_SEP = "/"
local FUNC_SEP = "_"


------------------------------------------------------------------------------

function compilefile(filename)
   local string = assert(loadfile(filename))
   fh = assert(io.open(filename,"wb"))
   fh:write(string.dump(string))
   fh:close()
end

function compilestring(str)
   return string.dump(assert(loadstring(str)))
end

------------------------------------------------------------------------------

local function getbytecode(filepath,path_sep)
   filepath = filepath:gsub("%?",path_sep)
   local file = io.open(filepath)
   assert(file,"could not open "..filepath)
   file:close()
   return string.dump(assert(loadfile(filepath)))
end

function precompile(inputs, libname, outdir, path_sep)
   path_sep = path_sep or PATH_SEP
   local outc = assert(io.open(outdir..path_sep..libname.."_lmod.c", "w"))
   local outh = assert(io.open(outdir..path_sep..libname.."_lmod.h", "w"))

   local guard = libname:upper():gsub("[^%w]", "_")

-- write <libname>.h file

   outh:write([[
#ifndef __]],guard,[[__
#define __]],guard,[[__

#include <lua.h>

#ifndef ]],PREFIX,[[ 
#define ]],PREFIX,[[ 
#endif

]])
   
-- write <libname>.c file

   outc:write([[
#include <lua.h>
#include <lauxlib.h>
#include "]],libname,[[_lmod.h"

]])

   local i = 0
   for input,v in pairs(inputs) do
      i = i + 1
      local bytecode = getbytecode(input,path_sep)
      outc:write("static const unsigned char B",i,"[]={\n")
      for j = 1, #bytecode do
	 outc:write(string.format("%3u,", bytecode:byte(j)))
	 if j % 20 == 0 then outc:write("\n") end
      end
      outc:write("\n};\n\n")
   end

   local i = 0
   for k,input in pairs(inputs) do
      i = i + 1
      local func = input:gsub("%.", FUNC_SEP)
      outh:write(PREFIX," int luaopen_",func,"(lua_State *L);\n")
      outc:write(PREFIX,[[ int luaopen_]],func,[[(lua_State *L) {
  int arg = lua_gettop(L);
  luaL_loadbuffer(L,(const char*)B]],i,[[,sizeof(B]],i,[[),"]],input,[[");
  lua_insert(L,1);
  lua_call(L,arg,1);
  return 1;
}

]])
   end

   outh:write([[

#endif /* __]],guard,[[__ */

]])

   outh:close()
   outc:close()

end

------------------------------------------------------------------------------

function preloader(inputs, libname, ldrname, outdir, path_sep)

local funcname = "luapreload_"..libname

-- write <ldrname>.h file

local outh = assert(io.open(outdir..path_sep..ldrname..".h", "w"))

local guard = libname:upper():gsub("[^%w]", "_")

outh:write([[
#ifndef __]],guard,[[__
#define __]],guard,[[__

#ifndef ]],PREFIX,[[ 
#define ]],PREFIX,[[ 
#endif

]],PREFIX,[[ int ]],funcname,[[(lua_State *L);

#endif /* __]],guard,[[__ */
]])
outh:close()


-- write <ldrname>.c file

local outc = assert(io.open(outdir..path_sep..ldrname..".c", "w"))
outc:write([[
#include <lua.h>
#include <lauxlib.h>

]])

for i, input in ipairs(inputs) do
   outc:write('int luaopen_',input:gsub("%.", FUNC_SEP),'(lua_State*);\n')
end

outc:write([[
#include "]],libname,[[_lmod.h"

]],PREFIX,[[ int ]],funcname,[[(lua_State *L) {
	luaL_findtable(L, LUA_GLOBALSINDEX, "package.preload", ]], #inputs, [[);
	
]])
local code = [[
	lua_pushcfunction(L, luaopen_%s);
	lua_setfield(L, -2, "%s");
]]
for i, input in ipairs(inputs) do
   outc:write(code:format(input:gsub("%.", FUNC_SEP), input))
end
outc:write([[
	
	lua_pop(L, 1);
	return 0;
}
]])

outc:close()

end

------------------------------------------------------------------------------


lib = "cfg"
list = { ["cfg.lua"] = "cfg" }

precompile(list,lib,".","/")

