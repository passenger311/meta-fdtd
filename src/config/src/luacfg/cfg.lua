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


local _G,print,pairs,ipairs = _G,print,pairs,ipairs
local geo = require "geo"

module("cfg")

function Material(args)
   print("NEW MATERIAL")
end














-- inject "geo" and "cfg" modules into the global namespace
for k,v in _G.pairs(_M) do _G[k] = v end
for k,v in _G.pairs(geo) do _G[k] = v end
