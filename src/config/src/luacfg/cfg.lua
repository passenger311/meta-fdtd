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


local _G,print,pairs,ipairs,type,assert,setmetatable,table,string,io,tostring,
   unpack = 
   _G,print,pairs,ipairs,type,assert,setmetatable,table,string,io,tostring,
   unpack
local geo = require "geo"
local math = math

module("cfg")

-- Config

ConfigMethods = {}
ConfigMethods.__index = ConfigMethods
 

---------------------------------------------------------------------------
-- BLOCK DEFINITIONS
---------------------------------------------------------------------------

local scenes, pscenes = {}, {}

function CONFIG(parms)
   local tab = {}
   tab.mat = {} -- subtables
   tab.src = {}
   tab.diag = {}
   tab.lumped = {}
   if parms.scenes == nil then parms.scenes = true end
   tab.scenes_on = parms.scenes
   return setmetatable(tab,ConfigMethods)
end

-- GRID Config-Block definition

function ConfigMethods:GRID(parms)
   self.grid = { block="GRID" }
   self.grid.irange = parms.irange or { 0, 0 };
   self.grid.jrange = parms.jrange or { 0, 0 };
   self.grid.krange = parms.krange or { 0, 0 };
   self.grid.dim =  parms.dim or 1;
   self.grid.partition = parms.partition or { 0,1 }
   if type(parms.ncyc) == 'table' then
      self.grid.ncyc = parms.ncyc
   else		       
      self.grid.ncyc = { 0, parms.ncyc or 100 }
   end	
   self.grid.dt = parms.dt or 0.9999;
   self.grid.dx = parms.dx or { 1., 1., 1., 1. }
end

-- CHECKPOINT Checkpoint-Block definition

function ConfigMethods:CHECKPOINT(parms)
   self.checkpoint = { block="CHECKPOINT" }
   self.checkpoint.load = parms.load or false;
   self.checkpoint.save = parms.save or false;
   self.checkpoint.detail = parms.detail or 3;
end

-- FDTD Config-Block definition

function ConfigMethods:FDTD(parms)
   self.fdtd = { block="FDTD" }
   for k,v in ipairs(parms) do self.fdtd[k] = v end
end

-- EPSILON Sub-Block definition

function EPSILON(parms)
   local EPSILON = { block="EPSILON" }
   for k,v in pairs(parms) do EPSILON[k] = v end 
   return EPSILON
end

-- BOUND Config-Block definition

function ConfigMethods:BOUND(parms)
   self.bound = { block = "BOUND" }
   for k,v in pairs(parms) do self.bound[k] = v end 
   self.bound.config = parms.config or { 0,0,0,0,0,0 }
   for i = 1,6 do 
      if not self.bound.config[i] then self.bound.config[i] = 0 end 
   end
end

-- PML Sub-Block definition

function PML(parms)
   local PML = { block = "PML" }
   for k,v in pairs(parms) do PML[k] = v end 
   PML.cells = parms.cells or 11
   PML.pot = parms.pot or 3.2
   PML.sigma = parms.sigma or 1.94444444444444
   PML.kappa = parms.kappa or 1.1
   return PML
end

-- REG Sub-Block definition

function REG(parms) 
   local REG = { block = "REG" }
   for k,v in pairs(parms) do REG[k] = v end 
   return REG
end

-- LOAD Sub-Block definition

function LOAD(parms) 
   local LOAD = { block = "LOAD" }
   for k,v in pairs(parms) do LOAD[k] = v end
   return LOAD
end


-- (NEW!) LOAD_GEO Sub-Block definition

function LOAD_GEO(parms) 
   local LOAD_GEO = { block = "LOAD_GEO" }
   for k,v in pairs(parms) do LOAD_GEO[k] = v end
   assert(scenes[parms[1]] ~=nil, "LOAD_GEO{} <name>="..parms[1].." does not exist!")
   LOAD_GEO.file = "geo_"..tostring(parms[1])..".in"
   return LOAD_GEO
end

-- (NEW!) CREATE_GEO Sub-Block definition

function ConfigMethods:CREATE_GEO(parms) 
   assert(scenes[parms[1]] == nil, "LOAD_GEO{} <name>="..parms[1].." is already in use!")
   local file = "geo_"..tostring(parms[1])..".in"
   scenes[parms[1]] = 1
   if not self.scenes_on or parms.on == false then return end 
   local geo_fh = geo.FileIN{file,comps=parms.comps};
   assert(parms.scene and parms.grid,"CREATE_GEO{} must define <grid> and <scene>")
   print("processing "..file)
   parms.grid:write{geo_fh, parms.scene, method=parms.method, silent=parms.silent}
   print(" done.\n")
end


-- (NEW!) CREATE_PREVIEW Sub-Block definition

function ConfigMethods:CREATE_PREVIEW(parms) 
   assert(pscenes[parms[1]] == nil, "SCENE{} <name>="..parms[1].." is already in use!")
   local file = "preview_"..tostring(parms[1])..".vtk"
   pscenes[parms[1]] = 1
   if not self.scenes_on or parms.on == false then return end 
   local geo_fh = geo.FileVTK{file,comps=parms.comps};
   assert(parms.scene and parms.grid,"CREATE_PREVIEW{} must define <grid> and <scene>")
   print("processing "..file)
   parms.grid:write{geo_fh, parms.scene, method=parms.method, silent=parms.silent}
   print(" done.\n")
end


-- POINT Sub-Block definition

function POINT(parms) 
   local POINT = { block = "POINT" }
   for k,v in pairs(parms) do POINT[k] = v end
   return POINT
end

-- BOX Sub-Block definition

function BOX(parms) 
   local BOX = { block = "BOX" }
   for k,v in pairs(parms) do BOX[k] = v end
   return BOX
end

-- OUT Sub-Block definition

function OUT(parms) 
   local OUT = { block = "OUT" }
   for k,v in pairs(parms) do OUT[k] = v end
   OUT.file = parms.file or { "VTK", "?" }
   OUT.type = parms.type or { "?", "N", ".T." }
   OUT.type[2] = OUT.type[2] or "N"
   OUT.type[3] = OUT.type[3] or ".T."
   OUT.time = parms.time or { 0,-1, 1 }  
   return OUT
end

-- MAT Config-Block definition

function ConfigMethods:MAT(parms)
   local MAT = { block = "MAT" }
   for k,v in pairs(parms) do MAT[k] = v end
   table.insert(self.mat, MAT)   
end

--(MAT)DRUDE Sub-Block definition

function DRUDE(parms) 
   local DRUDE = { block = "LORENTZ" } -- DRUDE implemented as a LORENTZ
   DRUDE.a = 1
   DRUDE.b = parms.gammapl or 0
   DRUDE.c = 0
   DRUDE.d = parms.invlambdapl^2
--   DRUDE.invlambdapl = parms.invlambdapl
--   DRUDE.gammapl = parms.gammapl or 0
--   DRUDE.order = parms.order or 1
   return DRUDE
end

-- (MAT)LHM Sub-Block definition

function LHM(parms) 
   local LHM = { block = "LHM" }
   LHM.invlambdapl = parms.invlambdapl
   LHM.gammapl = parms.gammapl or 0
   LHM.order = parms.order or 1
   return LHM
end

-- (MAT)LHM Sub-Block definition

function LHMGRAD(parms) 
   local LHMGRAD = { block = "LHMGRAD" }
   LHMGRAD.file = parms.file
   LHMGRAD.from = parms.from or { 0, 0, 0 }
   LHMGRAD.to = parms.to or { 0, 0, 0 }
   return LHMGRAD
end

-- (MAT)PEC Sub-Block definition

function PEC(parms) 
   local PEC = { block = "PEC" }
   return PEC
end

-- (MAT)LORENTZ Sub-Block definition

function LORENTZ(parms) 
   local LORENTZ = { block = "LORENTZ" } -- more generic LORENTZ
   LORENTZ.a = 1
   LORENTZ.b = 2 * ( parms.gammal or 0 )
   LORENTZ.c = parms.invlambdal^2
   LORENTZ.d = parms.invlambdal^2 * ( parms.deltaepsl or 1e-3 )
--   LORENTZ.invlambdal = parms.invlambdal
--   LORENTZ.gammal = parms.gammal or 0
--   LORENTZ.deltaepsl = parms.deltaepsl or 1e-3
   return LORENTZ
end

-- (MAT)DEBYE Sub-Block definition

function DEBYE(parms) 
   local DEBYE = { block = "LORENTZ" } -- DEBYE inplemented as LORENTZ
   DEBYE.a = 0
   DEBYE.b = (2*math.pi)^2 * (parms.taud or 0.) --3.14159265358979323846
   DEBYE.c = 1
   DEBYE.d = parms.deltaepsd or 1e-3
--   DEBYE.taud = parms.taud or 0.
--   DEBYE.deltaepsd = parms.deltaepsd or 1e-3
   return DEBYE
end

-- (MAT)CHI3 Sub-Block definition

function CHI3(parms) 
   local CHI3 = { block = "CHI3" }
   CHI3.invlambdar = parms.invlambdar
   CHI3.gammar = parms.gammar
   CHI3.chi3r = parms.chi3r
   CHI3.chi3k = parms.chi3k or 0
   CHI3.maxit = parms.maxit or 0
   return CHI3
end


-- (MAT)BLOCH Sub-Block definition

function BLOCH(parms) 
   local BLOCH = { block = "BLOCH" }
   BLOCH.invlambdal = parms.invlambdal
   BLOCH.gammal = parms.gammal or 0
   BLOCH.dipole = parms.dipole or { 0.1, 0.1, 0.1 }
   BLOCH.carrier = parms.carrier or { 1., 0.5 }
   BLOCH.gammanr = parms.gammanr or 0
   BLOCH.pump = parms.pump or 0
   return BLOCH
end

-- (MAT)FOURLVL Sub-Block definition

function FOURLVL(parms)
   local FOURLVL = { block = "FOURLVL" }
   FOURLVL.invlambdal = parms.invlambdal
   FOURLVL.gammal = parms.gammal or {0,0}
   FOURLVL.dipole12 = parms.dipole12 or 0.1
   FOURLVL.dipole03 = parms.dipole03 or 0.1
   FOURLVL.dens = parms.dens or 1
   FOURLVL.start = parms.start or {0,0,0}
   FOURLVL.gamma = parms.gamma or {0,0,0,0}
   return FOURLVL
end

-- (Mat) QW Sub-Block definition

function QW(parms)
   local QW = { block = "QW" }
   QW.ks = parms.ks or { 0., 1. }
   QW.kCount = parms.kCount or 100
   QW.bandgap = parms.bandgap or 1.5
   QW.masses = parms.masses or { 0.06, 0.3 }
   QW.dcv = parms.dcv or { 0.1 }
   QW.nr = parms.nr or 1
   QW.ninitial = parms.ninitial or 1.
   QW.T = parms.T or 300
   QW.pump = parms.pump or 0.0001
   QW.intragammas = parms.intragammas or { 0.001, 0.001, 0.001 }
   QW.gammamac = parms.gammamac or { 0.0001, 1.e-8, 1.e-8 }
   return QW
end

-- (MAT)THREELVL Sub-Block definition

function THREELVL(parms) 
   local THREELVL = { block = "THREELVL" }
   THREELVL.invlambda = parms.invlambda or {}
   THREELVL.gamma = parms.gamma or { 0., 0., 0. }
   THREELVL.sigma = parms.sigma or { 0., 0., 0. }
   THREELVL.mx = parms.mx or { {0,0}, {0,0}, {0,0} }
   THREELVL.my = parms.my or { {0,0}, {0,0}, {0,0} }
   THREELVL.mz = parms.mz or { {0,0}, {0,0}, {0,0} }
   THREELVL.densities = parms.densities or { 1., 0, 0 }
   THREELVL.n = parms.n or 1
   THREELVL.LFE = parms.LFE or 1
   THREELVL.epsLFE = parms.epsLFE or 1
   return THREELVL
end

-- (MAT)RANDTHREELVL Sub-Block definition

function RANDTHREELVL(parms)
   local RANDTHREELVL = { block = "RANDTHREELVL" }
   RANDTHREELVL.invlambda = parms.invlambda or {}
   RANDTHREELVL.gamma = parms.gamma or { 0., 0., 0. }
   RANDTHREELVL.sigma = parms.sigma or { 0., 0., 0. }
   RANDTHREELVL.Mvals = parms.Mvals or { 0.1, 0.1, 0.}
   RANDTHREELVL.angle = parms.angle or 0.
   RANDTHREELVL.densities = parms.densities or { 1., 0., 0. }
   RANDTHREELVL.n = parms.n or 1
   RANDTHREELVL.LFE = parms.LFE or 1
   RANDTHREELVL.seed = parms.seed or 1734920
   return RANDTHREELVL
end

-- (MAT)TWOLVL Sub-Block definition

function TWOLVL(parms) 
   local TWOLVL = { block = "TWOLVL" }
   TWOLVL.invlambda = parms.invlambda or {}
   TWOLVL.gamma = parms.gamma or 0.
   TWOLVL.sigma = parms.sigma or 0.
   TWOLVL.m = parms.m or { {0,0}, {0,0}, {0,0} }
   TWOLVL.inversion = parms.inversion or -1.
   TWOLVL.n = parms.n or 1
   TWOLVL.LFE = parms.LFE or 1
   TWOLVL.Lorentz = parms.Lorentz or 1
   return TWOLVL
end


-- SRC Config-Block definition

function ConfigMethods:SRC(parms)
   local SRC = { block = "SRC" }
   for k,v in pairs(parms) do SRC[k] = v end
   table.insert(self.src, SRC)   
end

--(SRC)HARDJ Sub-Block definition

function HARDJ(parms) 
   local HARDJ = { block = "HARDJ" }
   HARDJ.invlambda = parms.invlambda
   HARDJ.amplitude = parms.amplitude or 1.
   HARDJ.pulse = parms.pulse or { shape="Gaussian", width=0, 
				  offset=0, attack=0, sustain=0, decay=0 }
   HARDJ.planewave = parms.planewave or { on=false, phi=0, theta=0, psi=0, nrefr=1 }
   HARDJ.alpha = parms.alpha or 0
   HARDJ.domega = parms.domega or 0
   if not HARDJ.planewave.on then HARDJ.planewave.on=true end
   return HARDJ
end

--(SRC)TFSFINJ Sub-Block definition

function TFSFINJ(parms) 
   local TFSFINJ = { block = "TFSFINJ" }
   TFSFINJ.invlambda = parms.invlambda
   TFSFINJ.amplitude = parms.amplitude or 1.
   TFSFINJ.pulse = parms.pulse or { shape="Gaussian", width=0, 
				    offset=0, attack=0, sustain=0, decay=0 }
   TFSFINJ.planewave = parms.planewave or { on=true, phi=0, theta=0, psi=0, nrefr=1 }
   TFSFINJ.alpha = parms.alpha or 0
   TFSFINJ.domega = parms.domega or 0
   return TFSFINJ
end

--(SRC)TFSFBOX Sub-Block definition

function TFSFBOX(parms) 
   local TFSFBOX = { block = "TFSFBOX" }
   TFSFBOX.invlambda = parms.invlambda
   TFSFBOX.amplitude = parms.amplitude or 1.
   TFSFBOX.pulse = parms.pulse or { shape="Gaussian", width=0, 
				    offset=0, attack=0, sustain=0, decay=0 }
   TFSFBOX.planewave = parms.planewave or { on=true, phi=0, theta=0, psi=0, nrefr=1 }
   TFSFBOX.config = parms.config or { 1,1,1,1,1,1 }
   TFSFBOX.alpha = parms.alpha or 0
   TFSFBOX.domega = parms.domega or 0
   return TFSFBOX
end

-- DIAG Config-Block definition

function ConfigMethods:DIAG(parms)
   local DIAG = { block = "DIAG" }
   for k,v in pairs(parms) do DIAG[k] = v end
   table.insert(self.diag, DIAG)   
end


-- (DIAG)PSPEC Sub-Block definition

function PSPEC(parms) 
   local PSPEC = { block = "PSPEC" }
   PSPEC.file = parms.file
   PSPEC.reffile = parms.reffile or ""
   PSPEC.mode = parms.mode or "S"
   PSPEC.phasewrap = parms.phasewrap or { 0, 0 } 
   PSPEC.time = parms.time or { 0,-1, 1 }  
   PSPEC.polarize = parms.polarize or { phi=0, theta=0, psi=0 }
   return PSPEC
end

-- (DIAG)MODE Sub-Block definition

function MODE(parms) 
   local MODE = { block = "MODE" }
   MODE.time = parms.time or { 0,-1, 1 }  
   MODE.mode = parms.mode or "F"
   MODE.file = parms.file
   MODE.outfile = parms.outfile
   return MODE
end

-- (DIAG)EBAL Sub-Block definition

function EBAL(parms) 
   local EBAL = { block = "EBAL" }
   EBAL.time = parms.time or { 0,-1, 1 }  
   return EBAL
end

-- LUMPED Config-Block definition

function ConfigMethods:LUMPED(parms)
   local LUMPED = { block = "LUMPED" }
   for k,v in pairs(parms) do LUMPED[k] = v end
   LUMPED.type = parms.type or "C" 
   table.insert(self.lumped, LUMPED)   
end

---------------------------------------------------------------------------
-- OUTPUT METHODS
---------------------------------------------------------------------------

-- GRID Config-Block write

local function writeGRID(fh,GRID)
   assert(GRID and GRID.block == "GRID", "Expected GRID{}")
   fh:write("(GRID\n");
   fh:write("  ",GRID.dim," \t! dim\n");
   fh:write("  ",GRID.partition[1], " ", GRID.partition[2]," \t! partition (i of n)\n");
   fh:write("  ",GRID.ncyc[1], " ", GRID.ncyc[2], " \t! ncyc (# of timesteps)\n");
   fh:write("  ",GRID.dt," \t! dt\n");
   fh:write("  ",GRID.irange[1], " ", GRID.irange[2]," \t! irange = (ibeg,iend)\n");
   fh:write("  ",GRID.jrange[1], " ", GRID.jrange[2]," \t! jrange = (jbeg,jend)\n");
   fh:write("  ",GRID.krange[1], " ", GRID.krange[2]," \t! krange = (kbeg,kend)\n");
   fh:write("  ",GRID.dx[1], " ", GRID.dx[2]," ", GRID.dx[3], " ", GRID.dx[4], " \t! conv fac + step sizes\n");
   fh:write(")GRID\n\n");
end

-- CHECKPOINT Config-Block write

local function writeCHECKPOINT(fh,CHECKPOINT)
   assert(CHECKPOINT and CHECKPOINT.block == "CHECKPOINT", "Expected CHECKPOINT{}")
   fh:write("(CHECKPOINT\n");
   if ( CHECKPOINT.load ) then CHECKPOINT.load = ".T." else CHECKPOINT.load = ".F." end
   if ( CHECKPOINT.save ) then CHECKPOINT.save = ".T." else CHECKPOINT.save = ".F." end
   fh:write("  ",CHECKPOINT.load," \t! load\n");
   fh:write("  ",CHECKPOINT.save," \t! save\n");
   fh:write("  ",CHECKPOINT.detail," \t! detail level (1=field, 2=material, 3=field+material)\n");
   fh:write(")CHECKPOINT\n\n");
end

-- BOX Sub-Block write

local function writeBOX(fh,BOX)
   assert(BOX and BOX.block == "BOX", "Expected BOX{}")
   if BOX.on == false then return end  
   fh:write("    (BOX\n")
   for _,box in ipairs(BOX) do
      fh:write("\t")
      for _,val in ipairs(box) do
	 fh:write(val," ")
      end
      fh:write("\n");
   end
   fh:write("    )BOX\n")
end

-- POINT Sub-Block write

local function writePOINT(fh,POINT)
   assert(POINT and POINT.block == "POINT", "Expected POINT{}")
   if POINT.on == false then return end  
   fh:write("    (POINT\n")
   for _,point in ipairs(POINT) do
      fh:write("\t")
      for _,val in ipairs(point) do
	 fh:write(val," ")
      end
      fh:write("\n");
   end
   fh:write("    )POINT\n")
end

-- LOAD Sub-Block write

local function writeLOAD(fh,LOAD)
   assert(LOAD and LOAD.block == "LOAD", "Expected LOAD{}")
   if LOAD.on == false then return end  
   fh:write("    (LOAD\n")
   fh:write("\t",LOAD[1],"\t! file to load\n")
   fh:write("    )LOAD\n")
end

-- LOAD_GEO Sub-Block write

local function writeLOAD_GEO(fh,LOAD_GEO)
   assert(LOAD_GEO and LOAD_GEO.block == "LOAD_GEO", "Expected LOAD_GEO{}")
   if LOAD_GEO.on == false then return end  
   fh:write("    (LOAD\n")
   fh:write("\t",LOAD_GEO.file,"\t! scene file to load\n")
   fh:write("    )LOAD\n")
end

-- REG Sub-Block write

local function writeREG(fh,REG)
   assert(REG and REG.block == "REG", "Expected REG{}")
   if REG.on == false then return end  
   fh:write("  (REG\n")
   for i,v in ipairs(REG) do
      if v.block == "BOX" then writeBOX(fh,v) end
      if v.block == "POINT" then writePOINT(fh,v) end
      if v.block == "LOAD" then writeLOAD(fh,v) end
      if v.block == "LOAD_GEO" then writeLOAD_GEO(fh,v) end
   end
   if REG.auto then fh:write("    AUTO\t! auto loop mode\n") end
   if REG.mask then fh:write("    MASK\t! mask loop mode\n") end
   if REG.list then fh:write("    LIST\t! list loop mode\n") end
   fh:write("  )REG\n")
end

-- EPSILON Sub-Block write

local function writeEPSILON(fh,EPSILON)
   assert(EPSILON and EPSILON.block == "EPSILON", "Expected EPSILON{}")
   if EPSILON.on == false then return end  
   fh:write("(EPSILON\n")
   assert(EPSILON[1] and EPSILON[1].block == "REG","EPSILON{} must define REG{}")
   writeREG(fh, EPSILON[1])
   fh:write(")EPSILON\n")
end

-- OUT Sub-Block write

local function writeOUT(fh,OUT)
   assert(OUT and OUT.block == "OUT", "Expected OUT{}")
   assert(OUT[1] and OUT[1].block == "REG","OUT{} must define REG{}")
   if OUT.on == false then return end  
   fh:write("(OUT\n")
   fh:write("  ",OUT.file[1]," ",OUT.file[2],"\t! file type and name\n"); 
   fh:write("  ",OUT.type[1]," ",OUT.type[2]," ",OUT.type[3],"\t! file type and name\n"); 
   fh:write("  ",OUT.time[1]," ",OUT.time[2]," ",OUT.time[3],"\t! time\n");
   writeREG(fh, OUT[1])
   fh:write(")OUT\n")
end

-- FDTD Config-Block write

local function writeFDTD(fh,FDTD)
   assert(FDTD and FDTD.block == "FDTD", "Expected FDTD{}")
   fh:write("(FDTD\n");
   for _, sub in ipairs(FDTD) do 
      if sub.block then
	 if sub.block == "OUT" then 
	    writeOUT(fh,sub);
	 end
	 if sub.block == "EPSILON" then
	    writeEPSILON(fh,sub);
	 end
      end
   end
   fh:write(")FDTD\n\n");
end

-- PML Sub-Block write

local function writePML(fh,PML)
   assert(PML and PML.block == "PML", "Expected PML{}")
   fh:write("(PML\n")
   fh:write(PML.cells,"\t! # of cells\n");
   fh:write(PML.pot,"\t! pot parameter\n");
   fh:write(PML.sigma,"\t! sigma parameter\n");
   fh:write(PML.kappa,"\t! kappa parameter\n");
   fh:write(")PML\n")
end

-- BOUND Config-Block write

local function writeBOUND(fh,BOUND)
   assert(BOUND and BOUND.block == "BOUND", "Expected BOUND{}")
   fh:write("(BOUND\n")
   fh:write(BOUND.config[1]," ",BOUND.config[2]," ",BOUND.config[3]," ",
	    BOUND.config[4]," ",BOUND.config[5]," ",BOUND.config[6],
	    "\t! boundary config\n"); 
   for i,v in ipairs(BOUND) do
      if v.block == "PML" then writePML(fh,v) end
   end
   fh:write(")BOUND\n\n")
end

-- (MAT) MAT Sub-Block write

local writemat = {
--[[
   DRUDE = function(fh,DRUDE)
	      fh:write(1,"\t! a = 1\n")
	      fh:write(DRUDE.invlambdapl^2,"\t! invlambdapl [2 pi c]\n")
	      fh:write(DRUDE.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write(DRUDE.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
--]]
   LHM = function(fh,LHM)
	      fh:write(LHM.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write(LHM.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write(LHM.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
   LHMGRAD = function(fh,LHMGRAD)
	      fh:write(LHMGRAD.file,"\t! file to load\n")
	      fh:write(LHMGRAD.from[1]," ",
		       LHMGRAD.from[2]," ",
		       LHMGRAD.from[3],"\t! offset point\n")
	      fh:write(LHMGRAD.to[1]," ",
		       LHMGRAD.to[2]," ",
		       LHMGRAD.to[3],"\t! size vector\n")
	   end,
   LORENTZ = function(fh,LORENTZ)
	      fh:write(LORENTZ.a,"\t! a\n")
	      fh:write(LORENTZ.b,"\t! b (damping) [1/dt]\n")
	      fh:write(LORENTZ.c,"\t! c (invlambda^2) [(2 pi c)**2]\n")
	      fh:write(LORENTZ.d,"\t! d (coupling*invlambda^2)[(2 pi c)**2]\n")
	   end,
--[[
   DEBYE = function(fh,DEBYE)
	      fh:write(DEBYE.taud,"\t! taud (damping) [dt]\n")
	      fh:write(DEBYE.deltaepsd,"\t! deltaepsd (coupling coeff.)\n")
	   end,
--]]
   PEC = function(fh,PEC)
	 end,
   BLOCH = function(fh,BLOCH)
	      fh:write(BLOCH.invlambdal,"\t! invlambdal [2 pi c]\n")
	      fh:write(BLOCH.gammal,"\t! gammal (damping) [1/dt]\n")
	      fh:write(BLOCH.dipole[1], " " , BLOCH.dipole[2]," ", BLOCH.dipole[3],"\t! dipole matrix elem. [1/dx]\n")
	      fh:write(BLOCH.carrier[1]," ", BLOCH.carrier[2],
		       "\t! carrier numbers transp./initial  [1/(dx^3)]\n")
	      fh:write(BLOCH.gammanr,"\t! non-rad. recomb. [1/dt]\n")
	      fh:write(BLOCH.pump,"\t! pump rate [1/dt]\n")		       
	   end,
   FOURLVL = function(fh,FOURLVL)
              fh:write(FOURLVL.invlambdal[1]," ", FOURLVL.invlambdal[2],"\t! invlambdal [2 pi c]\n")
	      fh:write(FOURLVL.gammal[1]," ", FOURLVL.gammal[2],"\t! gammal [1/dt]\n")
	      fh:write(FOURLVL.dipole12,"\t! DPME 1<->2 [1/dx]\n")
	      fh:write(FOURLVL.dipole03,"\t! DPME 0<->3 [1/dx]\n")
	      fh:write(FOURLVL.dens,"\t! density of four lvl systems [1/dx^3]\n")
	      fh:write(FOURLVL.start[1]," ", FOURLVL.start[2]," ",FOURLVL.start[3],"\t! values for N3_0,N2_0,N1_0 [1/dx^3]\n")
	      fh:write(FOURLVL.gamma[1]," ",FOURLVL.gamma[2]," ",FOURLVL.gamma[3]," ",FOURLVL.gamma[4],"\t! non rad. recomb. rates [1/dt]\n")
	  end,
   QW = function(fh,QW)
              fh:write(QW.ks[1]," ", QW.ks[2],"\t! { k0, kmax } [1/dx]\n")
	      fh:write(QW.kCount,"\t! number of momentum states\n")
	      fh:write(QW.bandgap,"\t! bandgap in eV\n")
	      fh:write(QW.masses[1]," ", QW.masses[2],"\t! { me, mh} [m_0]\n")
	      fh:write(QW.dcv,"\t! dipole matrix element [dx]\n")
	      fh:write(QW.nr,"\t! number of quantum wells per dx\n")
	      fh:write(QW.ninitial,"\t! initial electron density [1/dx^2]\n")
	      fh:write(QW.T,"\t! temperature in Kelvin\n")
	      fh:write(QW.pump,"\t! pump rate [1/dt]\n")
	      fh:write(QW.intragammas[1]," ",QW.intragammas[2]," ",QW.intragammas[3],"\t! {gamma_p, gamma_e, gamma_h} [1/dt]\n")
	      fh:write(QW.gammamac[1]," ",QW.gammamac[2]," ",QW.gammamac[3],"\t! {gamma_nr, gamma_sp, gamma_aug} [1/dt]\n")
	  end,
   CHI3 = function(fh,CHI3)
	      fh:write(CHI3.invlambdar,"\t! invlambdar [2 pi c]\n")
	      fh:write(CHI3.gammar,"\t! gammar (damping) [1/dt]\n")
	      fh:write(CHI3.chi3r,"\t! chi3 for raman\n")
	      fh:write(CHI3.chi3k,"\t! chi3 for kerr\n")
	      fh:write(CHI3.maxit,"\t! maximum iterations for kerr\n")
	  end,
   THREELVL = function(fh,THREELVL)
	      fh:write(THREELVL.invlambda[1]," ",THREELVL.invlambda[2]," ", THREELVL.invlambda[3], "\n")
	      fh:write(THREELVL.gamma[1]," ",THREELVL.gamma[2]," ", THREELVL.gamma[3], "\n")
	      fh:write(THREELVL.sigma[1]," ",THREELVL.sigma[2]," ", THREELVL.sigma[3], "\n")
	      fh:write(
		       "(",THREELVL.mx[1][1],",",THREELVL.mx[1][2],") ",
		       "(",THREELVL.mx[2][1],",",THREELVL.mx[2][2],") ",
		       "(",THREELVL.mx[3][1],",",THREELVL.mx[3][2],") ",
		       "\t! mx dipole length [dx]\n")
	      fh:write(
		       "(",THREELVL.my[1][1],",",THREELVL.my[1][2],") ",
		       "(",THREELVL.my[2][1],",",THREELVL.my[2][2],") ",
		       "(",THREELVL.my[3][1],",",THREELVL.my[3][2],") ",
		       "\t! my dipole length [dx]\n")
	      fh:write(
		       "(",THREELVL.mz[1][1],",",THREELVL.mz[1][2],") ",
		       "(",THREELVL.mz[2][1],",",THREELVL.mz[2][2],") ",
		       "(",THREELVL.mz[3][1],",",THREELVL.mz[3][2],") ",
		       "\t! mz dipole length [dx]\n")
	      fh:write(THREELVL.densities[1]," ", THREELVL.densities[2]," ", THREELVL.densities[3], "\t! occup. densities []\n")
	      fh:write(THREELVL.n,"\t! systems per cell []\n")
	      fh:write(THREELVL.LFE,"\t! local field effect included?\n")
	      fh:write(THREELVL.epsLFE,"\t! local field effect due to epsilon included?\n")
	   end,
   RANDTHREELVL = function(fh,RANDTHREELVL)
              fh:write(RANDTHREELVL.invlambda[1]," ",RANDTHREELVL.invlambda[2]," ",RANDTHREELVL.invlambda[3], "\n")
	      fh:write(RANDTHREELVL.gamma[1]," ",RANDTHREELVL.gamma[2]," ",RANDTHREELVL.gamma[3], "\n")
	      fh:write(RANDTHREELVL.sigma[1]," ",RANDTHREELVL.sigma[2]," ",RANDTHREELVL.sigma[3], "\n")
	      fh:write(RANDTHREELVL.Mvals[1]," ",RANDTHREELVL.Mvals[2]," ",RANDTHREELVL.Mvals[3], "\n")
	      fh:write(RANDTHREELVL.angle, "\n")
	      fh:write(RANDTHREELVL.densities[1]," ",RANDTHREELVL.densities[2]," ",RANDTHREELVL.densities[3], "\n")
	      fh:write(RANDTHREELVL.n, "\n")
	      fh:write(RANDTHREELVL.LFE, "\n")
	      fh:write(RANDTHREELVL.seed, "\n")
	   end,
	  TWOLVL = function(fh,TWOLVL)
	      fh:write(TWOLVL.invlambda, "\n")
	      fh:write(TWOLVL.gamma, "\n")
	      fh:write(TWOLVL.sigma, "\n")
	      fh:write(
		       "(",TWOLVL.m[1][1],",",TWOLVL.m[1][2],") ",
		       "(",TWOLVL.m[2][1],",",TWOLVL.m[2][2],") ",
		       "(",TWOLVL.m[3][1],",",TWOLVL.m[3][2],") ",
		       "\t! m dipole length [dx]\n")
	      fh:write(TWOLVL.inversion, "\t! occupation inversion []\n")
	      fh:write(TWOLVL.n,"\t! systems per cell []\n")
	      fh:write(TWOLVL.LFE,"\t! local field effect included?\n")
		  fh:write(TWOLVL.Lorentz,"\t! Lorentz local field correction included?\n")
	   end,
}

-- MAT Config-Block write

local function writeMAT(fh,MAT)
   assert(MAT and MAT.block == "MAT", "Expected MAT{}")
   if MAT.on == false then return end  
   assert( MAT[1] and MAT[2] and MAT[1].block and MAT[2].block, "Bad MAT{} structure") 
   local type = string.upper(MAT[1].block)
   fh:write("(MAT"..type.."\n")
   assert(writemat[type], "MAT{} type "..type.." does not exist")
   writemat[type](fh,MAT[1])
   writeREG(fh, MAT[2])
   for _, sub in ipairs(MAT) do 
      if sub.block then
	 if sub.block == "OUT" then 
	    writeOUT(fh,sub);
	 end
      end
   end  
   fh:write(")MAT"..type.."\n\n");
end

-- (SRC) SRC Sub-Block write

local writesrc = {
   HARDJ = function(fh,HARDJ)
	      fh:write(HARDJ.invlambda," \t! invlambda [2 pi c]\n")
	      fh:write(HARDJ.amplitude," \t! amplitude []\n")
	      fh:write(HARDJ.pulse.shape," \t! pulse shape\n")
	      fh:write(HARDJ.pulse.width," \t! pulse width [hwhm]\n")
	      fh:write(HARDJ.pulse.offset or 0," ",
		       HARDJ.pulse.attack or 0," ",
		       HARDJ.pulse.sustain or 0," ",
		       HARDJ.pulse.decay or 0,"\t! offset attack sustain decay [dt]\n")
	      if HARDJ.planewave.on then
		 fh:write(".T."," \t! plane wave mode?\n ")
	      else
		 fh:write(".F."," \t! plane wave mode?\n ")
	      end
	      fh:write(HARDJ.planewave.phi, " ",
		       HARDJ.planewave.theta, " ",
		       HARDJ.planewave.psi, " ",
		       HARDJ.planewave.nrefr, " \t! planewave: phi, theta, psi, nrefr\n")
	      fh:write(HARDJ.alpha," \t! carrier envelope phase \n")
	      fh:write(HARDJ.domega, " \t! linear chirp \n")
	   end,
   TFSFINJ = function(fh,TFSFINJ)
	      fh:write(TFSFINJ.invlambda," \t! invlambda [2 pi c]\n")
	      fh:write(TFSFINJ.amplitude," \t! amplitude []\n")
	      fh:write(TFSFINJ.pulse.shape," \t! pulse shape\n")
	      fh:write(TFSFINJ.pulse.width," \t! pulse width [hwhm]\n")
	      fh:write(TFSFINJ.pulse.offset or 0," ",
		       TFSFINJ.pulse.attack or 0," ",
		       TFSFINJ.pulse.sustain or 0," ",
		       TFSFINJ.pulse.decay or 0," \t! offset attack sustain decay [dt]\n")
	      fh:write(TFSFINJ.planewave.phi, " ",
		       TFSFINJ.planewave.theta, " ",
		       TFSFINJ.planewave.psi, " ",
		       TFSFINJ.planewave.nrefr, " \t! planewave: phi, theta, psi, nrefr\n")
	      fh:write(TFSFINJ.alpha," \t! relative phase between envelope and carrier \n")
	      fh:write(TFSFINJ.domega, " \t! linear chirp \n")
	   end,
   TFSFBOX = function(fh,TFSFBOX)
	      fh:write(TFSFBOX.invlambda," \t! invlambda [2 pi c]\n")
	      fh:write(TFSFBOX.amplitude," \t! amplitude []\n")
	      fh:write(TFSFBOX.pulse.shape," \t! pulse shape\n")
	      fh:write(TFSFBOX.pulse.width," \t! pulse width [hwhm]\n")
	      fh:write(TFSFBOX.pulse.offset or 0," ",
		       TFSFBOX.pulse.attack or 0," ",
		       TFSFBOX.pulse.sustain or 0," ",
		       TFSFBOX.pulse.decay or 0," \t! offset attack sustain decay [dt]\n")
	      fh:write(TFSFBOX.planewave.phi, " ",
		       TFSFBOX.planewave.theta, " ",
		       TFSFBOX.planewave.psi, " ",
		       TFSFBOX.planewave.nrefr, " \t! planewave: phi, theta, psi, nrefr\n")
	      fh:write(TFSFBOX.config[1]," ", TFSFBOX.config[2]," ",
		        TFSFBOX.config[3]," ", TFSFBOX.config[4]," ",
			TFSFBOX.config[5]," ", TFSFBOX.config[6]," \t! active planes\n")
              fh:write(TFSFBOX.alpha," \t! carrier envelope phase \n")
	      fh:write(TFSFBOX.domega, " \t! linear chirp \n")
	   end,

}

-- SRC Config-Block write

local function writeSRC(fh,SRC)
   assert(SRC and SRC.block == "SRC", "Expected SRC{}")
   if SRC.on == false then return end  
   assert( SRC[1] and SRC[2] and SRC[1].block and SRC[2].block, "Bad SRC{} structure") 
   local type = string.upper(SRC[1].block)
   fh:write("(SRC"..type.."\n")
   assert(writesrc[type], "SRC{} type "..type.." does not exist")
   writesrc[type](fh,SRC[1])
   writeREG(fh, SRC[2])
   for _, sub in ipairs(SRC) do 
      if sub.block then
	 if sub.block == "OUT" then 
	    writeOUT(fh,sub);
	 end
      end
   end  
   fh:write(")SRC"..type.."\n\n");
end

-- (DIAG) DIAG Sub-Block write

local writediag = {
   PSPEC = function(fh,PSPEC)
	      fh:write(PSPEC.file," \t! filename (.pspec)\n")
	      fh:write(PSPEC.mode, " ", PSPEC.reffile," \t! mode ( S,Ecs,Hcs,Eap,Hap) and [ref. file]\n")
	      fh:write(PSPEC.phasewrap[1]," ", PSPEC.phasewrap[2], " \t! unwrap phase forward / backward\n")
	      fh:write(PSPEC.time[1]," ",PSPEC.time[2]," ",PSPEC.time[3]," \t! time window [from to step]\n")
	      fh:write(PSPEC.polarize.phi, " ",
		       PSPEC.polarize.theta, " ",
		       PSPEC.polarize.psi, " \t! polarize: phi, theta, psi\n")
	   end,
   MODE =  function(fh,MODE)
	      fh:write(MODE.file," \t! filename, frequencies\n")
	      fh:write(MODE.outfile," \t! output filename\n")
	      fh:write(MODE.mode, " \t! mode ( F, En )\n")
	      fh:write(MODE.time[1]," ",MODE.time[2]," ",MODE.time[3]," \t! time window [from to step]\n")

	   end,
   EBAL = function(fh,EBAL)
	     fh:write(EBAL.time[1]," ",EBAL.time[2]," ",EBAL.time[3]," \t! time window [from to step]\n")
	  end
}


-- DIAG  Config-Block write

local function writeDIAG(fh,DIAG)
   assert(DIAG and DIAG.block == "DIAG", "Expected DIAG{}")
   if DIAG.on == false then return end  
   assert( DIAG[1] and DIAG[2] and DIAG[1].block and DIAG[2].block, "Bad DIAG{} structure") 
   local type = string.upper(DIAG[1].block)
   fh:write("(DIAG"..type.."\n")
   assert(writediag[type], "DIAG{} type "..type.." does not exist")
   writediag[type](fh,DIAG[1])
   writeREG(fh, DIAG[2])
   for _, sub in ipairs(DIAG) do 
      if sub.block then
	 if sub.block == "OUT" then 
	    writeOUT(fh,sub);
	 end
      end
   end  
   fh:write(")DIAG"..type.."\n\n");
end


-- LUMPED Config-Block write

local function writeLUMPED(fh,LUMPED)
   assert(LUMPED and LUMPED.block == "LUMPED", "Expected LUMPED{}")
   if LUMPED.on == false then return end  
   assert( LUMPED[1] , "Bad LUMPED{} structure") 
   local type = string.upper(LUMPED.type)
   fh:write("(LUMPED\n")
   fh:write(type.."\n");
   writeREG(fh, LUMPED[1])
   fh:write(")LUMPED\n\n");
end

-- CREATE configuration file

function ConfigMethods:CREATE()
   local part = self.grid.partition[1];
   local filename = "config."..tostring(part)..".in"
   local fh = io.open(filename,"w");
   fh:write("\n! ------ BEGIN [",filename,"] file generated by <luacfg>\n\n");
   writeGRID(fh,self.grid);
   if ( self.checkpoint ) then writeCHECKPOINT(fh,self.checkpoint) end
   writeFDTD(fh,self.fdtd);
   writeBOUND(fh,self.bound);
   for _, src in ipairs(self.src) do writeSRC(fh,src) end
   for _, mat in ipairs(self.mat) do writeMAT(fh,mat) end
   for _, diag in ipairs(self.diag) do writeDIAG(fh,diag) end
   for _, lumped in ipairs(self.lumped) do writeLUMPED(fh,lumped) end
    fh:write("! ------ END [",filename,"] \n\n");
   fh:close();
end


---------------------------------------------------------------------------
-- FINISH
---------------------------------------------------------------------------

-- move things to global namespace

for k,v in pairs(_M) do _G[k] = v end
for k,v in pairs(geo) do _G[k] = v end
