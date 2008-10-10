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
local geo_core = require "geo_core"

module("geo")

-- Object 

CObjectMethods = {}
CObjectMethods.__index = CObjectMethods

function CObjectMethods:__gc()
   geo_core.CObject_destroy(self)
end

function Sphere(parms)
   local obj = geo_core.Sphere_create(parms) 
   setmetatable(obj, CObjectMethods)
   return obj
end
function BezierPrism(parms)
   local obj = CBezierPrism_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function Box(parms)
   local obj = CBox_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function ConvexPrism(parms)
   local obj = CConvexPrism_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function SimplePrism(parms)
   local obj = CSimplePrism_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function Cylinder(parms)
   local obj = CCylinder_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function Ellipsoid(parms)
   local obj = CEllipsoid_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function BinaryAndNot(parms)
   local obj = CLogicAndNotObject_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function BinaryAnd(parms)
   local obj = CLogicAndObject_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function BinaryOr(parms)
   local obj = CLogicOrObject_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function BinaryXOr(parms)
   local obj = CLogicXOrObject_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj
end
function OrContainer(parms)
   local obj = CLogicOrContainer_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj   
end
function UnaryRotationZ(parms) 
   local obj = CSimpleRotationZ_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj   
end
function UnaryTransform(parms)
   local obj = CSimpleTransform_create(parms)
   setmetatable(obj, CObjectMethods)
   return obj   
end

-- Scene 

SceneMethods = {}
SceneMethods.__index = SceneMethods

function SceneMethods:__gc()
   geo_core.Scene_destroy(self)
end

function Scene(...)
   local scene = geo_core.Scene_create() 
   setmetatable(scene, SceneMethods)
   for i = 1, select("#",...) do 
      local obj = select(i,...) 
      assert(getmetatable(obj) == ObjectMethods)
      scene:add(obj)
   end
   return scene
end

-- Grid

GridMethods = {}
GridMethods.__index = GridMethods

function GridMethods:__gc()
   geo_core.Grid_destroy(self)
end

function Grid(scene, parms)
   assert(getmetatable(scene) == SceneMethods)
   local grid = geo_core.Grid_create(scene,parms) 
   setmetatable(grid, GridMethods)
   return grid
end
