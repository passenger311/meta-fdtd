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
 

----------------- definitions


function Config()
   local tab = {}
   tab.MAT = {} -- subtables
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
   self.GRID = { block="GRID" }
   self.GRID.irange = parms.irange or { 0, 0 };
   self.GRID.jrange = parms.jrange or { 0, 0 };
   self.GRID.krange = parms.krange or { 0, 0 };
   self.GRID.dim =  parms.dim or 1;
   self.GRID.partition = parms.partition or { 1,1 };
   self.GRID.ncyc = parms.ncyc or 100;
   self.GRID.dt = parms.ncyc or 0.9999;
end

-- Fdtd Block definition

function ConfigMethods:Fdtd(parms)
   self.FDTD = { block="FDTD" }
   for i,v in ipairs(parms) do self.FDTD[i] = parms[i] end
end

-- Epsilon Sub-Block definition

function Epsilon(parms)
   local EPSILON = { block="EPSILON" }
   for i,v in ipairs(parms) do EPSILON[i] = parms[i] end 
   return EPSILON
end

-- Bound Block definition

function ConfigMethods:Bound(parms)
   self.BOUND = { block = "BOUND" }
   self.BOUND.config = parms.config or { 0,0,0,0,0,0 }
   for i,v in ipairs(parms) do self.BOUND[i] = parms[i] end
end

-- Pml Sub-Block definition

function Pml(parms)
   local PML = { block = "PML" }
   PML.cells = parms.cells or 11
   PML.pot = parms.pot or 3.2
   PML.sigma = parms.sigma or 1.94444444444444
   PML.kappa = parms.kappa or 1.1
   for i,v in ipairs(parms) do PML[i] = parms[i] end
   return PML
end

-- Reg Sub-Block definition

function Reg(parms) 
   local REG = { block = "REG" }
   REG.mask = parms.mask
   REG.list = parms.list
   REG.auto = parms.auto
   for i,v in ipairs(parms) do REG[i] = parms[i] end
   return REG
end

-- Load Sub-Block definition

function Load(parms) 
   local LOAD = { block = "LOAD" }
   for i,v in ipairs(parms) do LOAD[i] = parms[i] end
   return LOAD
end

-- Point Sub-Block definition

function Point(parms) 
   local POINT = { block = "POINT" }
   for i,v in ipairs(parms) do POINT[i] = parms[i] end
   return POINT
end

-- Box Sub-Block definition

function Box(parms) 
   local BOX = { block = "BOX" }
   for i,v in ipairs(parms) do BOX[i] = parms[i] end
   return BOX
end

-- Out Sub-Block definition

function Out(parms) 
   local OUT = { block = "OUT" }
   OUT.file = parms.file or { "VTK", "undef" }
   OUT.type = parms.type or { "undef", "undef" }
   OUT.time = parms.time or { 0,-1, 1 }  
   for i,v in ipairs(parms) do OUT[i] = parms[i] end
   return OUT
end

-- Mat Block definition

function ConfigMethods:Mat(parms)
   local MAT = { block = "MAT" }
   for i,v in ipairs(parms) do MAT[i] = parms[i] end
   table.insert(self.MAT, MAT)   
end

-- Drude Sub-Block definition

function Drude(parms) 
   local DRUDE = { block = "DRUDE" }
   DRUDE.invlambdapl = parms.invlambdapl
   DRUDE.gammapl = parms.gammapl or 0
   DRUDE.order = parms.order or 1
   return DRUDE
end


-- Src Block definition

function ConfigMethods:Src(parms)
   local SRC = { block = "SRC" }
   for i,v in ipairs(parms) do SRC[i] = parms[i] end
   table.insert(self.SRC, SRC)   
end




----------------- end of definitions; processing follows


function ConfigMethods:write()
   

end







----------------- move things into global namespace

for k,v in pairs(_M) do _G[k] = v end