------------------------------------------------------------------------------
-- /file: cfg.lua
--
-- /author: J.Hamm 8/10/2008
--
-- /description:
--
-- The "cfg" module is the lua part of the configuration environment.
-- It includes the native C "geo" module and injects all entries in the
-- "geo" and "cfg" module tables into the global namespace.
--
------------------------------------------------------------------------------


local _G,print,pairs,ipairs,type,assert,setmetatable,table = 
   _G,print,pairs,ipairs,type,assert,setmetatable,table
local geo = require "geo"

module("cfg")

-- Config

ConfigMethods = {}
ConfigMethods.__index = ConfigMethods
 
function ConfigMethods:write(fh)
   -- TODO: write configuration to file

end

function Config()
   return setmetatable({},ConfigMethods)
end

-- write configuration to write handle

function CONFIG:write(fh)
   assert(self.GRID and self.FDTD)
   self.GRID:write(fh)
   self.FDTD:write(fh)
   self.BOUND:write(fh)
   self.SRC:write(fh)
   self.MAT:write(fh)
   self.DIAG:write(fh)
end

-- create a GRID block within configuration

GRID = {} 
GRID.__index = GRID

function GRID:write(fh)
   fh:write("(GRID",[[
]], self.is, self.ie, [[
]], self.js, self.je, [[
]], self.ks, self.ke, [[
]],")GRID")
end

function CONFIG:GRID(arg)   
   assert(arg.dim == 1 or arg.dim == 2 or arg.dim == 3)
   assert(type(arg.is) == "number" and type(arg.ie) == "number")
   assert(type(arg.js) == "number" and type(arg.je) == "number")
   assert(type(arg.ks) == "number" and type(arg.ke) == "number")
   grid = setmetatable({},GRID)
   table.insert(self.grids,grid)
   grid.dim = arg.dim
   grid.is,grid.ie = arg.is,arg.ie
   grid.js,grid.je = arg.js,arg.je
   grid.ks,grid.ke = arg.ks,arg.ke
   self.GRID = grid
end



-- create a FDTD block within configuration



-- Geometry objects:

function CONFIG:Scene(arg)
   scene = geo.Scene_create() 
   table.insert(self.scenes,scene)
   return scene
end

function CONFIG:Grid(arg)
   grid = geo.Grid_create() 
   table.insert(self.grids,grid)
   return grid
end
