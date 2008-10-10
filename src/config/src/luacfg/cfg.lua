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

function ConfigMethods:write(fh)
 
end

