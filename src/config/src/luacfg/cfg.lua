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


local _G,print,pairs,ipairs,type,assert,setmetatable,table,string = 
   _G,print,pairs,ipairs,type,assert,setmetatable,table,string
local geo = require "geo"

module("cfg")

-- Config

ConfigMethods = {}
ConfigMethods.__index = ConfigMethods
 
function ConfigMethods:write(fh)
   -- TODO: write configuration to file
end

function Config()
   local tab = {}
   tab.MAT = {}
   tab.SRC = {}
   tab.DIAG = {}
   return setmetatable(tab,ConfigMethods)
end

-- Reg Table

function Reg(parms)
   local reg = {}

   return reg
end

-- Out Table 

function Out(parms)
   local out = {}

   return out
end

-- Grid Block definition

function ConfigMethods:Grid(parms)
   self.GRID = {}
   
   
end

-- Fdtd Block definition

function ConfigMethods:Fdtd(parms)
   self.FDTD = {}
   
   
end

-- Bound Block definition

function ConfigMethods:Fdtd(parms)
   self.BOUND = {}
   
   
end

-- Mat Block definition

function ConfigMethods:Mat(name, parms)
   assert(type(name) == "string")
   name = string.toupper(name)
   local mat = {}
   table.insert(self.MAT,mat)

end

-- Src Block definition

function ConfigMethods:Src(name, parms)
   assert(type(name) == "string")
   name = string.toupper(name)
   local src = {}
   table.insert(self.MAT,src)
   

end
