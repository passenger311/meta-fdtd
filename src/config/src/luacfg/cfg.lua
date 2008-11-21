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
   self.grid.partition = parms.partition or { 0,1 };
   self.grid.ncyc = parms.ncyc or 100;
   self.grid.dt = parms.dt or 0.9999;
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
   local DRUDE = { block = "DRUDE" }
   DRUDE.invlambdapl = parms.invlambdapl
   DRUDE.gammapl = parms.gammapl or 0
   DRUDE.order = parms.order or 1
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
   LHMGRAD.point = parms.offset or { 0, 0, 0 }
   LHMGRAD.size = parms.size or { 0, 0, 0 }
   return LHMGRAD
end

-- (MAT)PEC Sub-Block definition

function PEC(parms) 
   local PEC = { block = "PEC" }
   return PEC
end

-- (MAT)LORENTZ Sub-Block definition

function LORENTZ(parms) 
   local LORENTZ = { block = "LORENTZ" }
   LORENTZ.invlambdal = parms.invlambdal
   LORENTZ.gammal = parms.gammal or 0
   LORENTZ.deltaepsl = parms.deltaepsl or 1e-3
   return LORENTZ
end

-- (MAT)DEBYE Sub-Block definition

function DEBYE(parms) 
   local DEBYE = { block = "DEBYE" }
   DEBYE.taud = parms.taud or 0.
   DEBYE.deltaepsd = parms.deltaepsd or 1e-3
   return DEBYE
end


-- (MAT)BLOCH Sub-Block definition

function BLOCH(parms) 
   local BLOCH = { block = "BLOCH" }
   BLOCH.invlambdal = parms.invlambdal
   BLOCH.gammal = parms.gammal or 0
   BLOCH.dipole = parms.dipole or { {1.,0.},{1.,0.},{1.,0.} }
   BLOCH.carrier = parms.carrier or { 1., 0.5 }
   BLOCH.gammanr = parms.gammanr or 0
   BLOCH.pump = parms.pump or 0
   BLOCH.satmodel = parms.satmodel or 0
   return BLOCH
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

-- (DIAG)EBAL Sub-Block definition

function EBAL(parms) 
   local EBAL = { block = "EBAL" }
   EBAL.time = parms.time or { 0,-1, 1 }  
   return EBAL
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
   fh:write("  ",GRID.ncyc," \t! ncyc (# of timesteps)\n");
   fh:write("  ",GRID.dt," \t! dt\n");
   fh:write("  ",GRID.irange[1], " ", GRID.irange[2]," \t! irange = (ibeg,iend)\n");
   fh:write("  ",GRID.jrange[1], " ", GRID.jrange[2]," \t! jrange = (jbeg,jend)\n");
   fh:write("  ",GRID.krange[1], " ", GRID.krange[2]," \t! krange = (kbeg,kend)\n");
   fh:write(")GRID\n\n");
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
   DRUDE = function(fh,DRUDE)
	      fh:write(DRUDE.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write(DRUDE.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write(DRUDE.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
   LHM = function(fh,LHM)
	      fh:write(LHM.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write(LHM.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write(LHM.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
   LHMGRAD = function(fh,LHMGRAD)
	      fh:write(LHMGRAD.file,"\t! file to load\n")
	      fh:write(LHMGRAD.point[1]," ",
		       LHMGRAD.point[2]," ",
		       LHMGRAD.point[3],"\t! offset point\n")
	      fh:write(LHMGRAD.size[1]," ",
		       LHMGRAD.size[2]," ",
		       LHMGRAD.size[3],"\t! size vector\n")
	   end,
   LORENTZ = function(fh,LORENTZ)
	      fh:write(LORENTZ.invlambdal,"\t! invlambdal [2 pi c]\n")
	      fh:write(LORENTZ.gammal,"\t! gammal (damping) [1/dt]\n")
	      fh:write(LORENTZ.deltaepsl,"\t! deltaepsl coupling coeff.\n")
	   end,
   DEBYE = function(fh,DEBYE)
	      fh:write(DEBYE.taud,"\t! taud [dt]\n")
	      fh:write(DEBYE.deltaepsd,"\t! deltaepsd coupling coeff.\n")
	   end,
   PEC = function(fh,PEC)
	 end,
   BLOCH = function(fh,BLOCH)
	      fh:write(BLOCH.invlambdal,"\t! invlambdal [2 pi c]\n")
	      fh:write(BLOCH.gammal,"\t! gammal (damping) [1/dt]\n")
	      fh:write(
		       "(",BLOCH.dipole[1][1],",",BLOCH.dipole[1][2],") ",
		       "(",BLOCH.dipole[2][1],",",BLOCH.dipole[2][2],") ",
		       "(",BLOCH.dipole[3][1],",",BLOCH.dipole[3][2],") ",
		       "\t! dipole matrix elem. []\n")
	      fh:write(BLOCH.carrier[1]," ", BLOCH.carrier[2],
		       "\t! carrier numbers transp./initial  []\n")
	      fh:write(BLOCH.gammanr,"\t! non-rad. recomb. [1/dt]\n")
	      fh:write(BLOCH.pump,"\t! pump rate [1/dt]\n")
	      fh:write(BLOCH.satmodel,"\t! sat.model 0=>(N-Ntr), 1=>Ntr*log(N/Ntr)\n")		       
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
		       HARDJ.pulse.decay or 0," \t! offset attack sustain decay [dt]\n")
	      if HARDJ.planewave.on then
		 fh:write(".T."," \t! plane wave mode?\n ")
	      else
		 fh:write(".F."," \t! plane wave mode?\n ")
	      end
	      fh:write(HARDJ.planewave.phi, " ",
		       HARDJ.planewave.theta, " ",
		       HARDJ.planewave.psi, " ",
		       HARDJ.planewave.nrefr, " \t! planewave: phi, theta, psi, nrefr\n")
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
   EBAL = function(fh,EBAL)
	     fh:write(EBAL.time[1]," ",PSPEC.time[2]," ",PSPEC.time[3]," \t! time window [from to step]\n")
	  end,
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

-- CREATE configuration file

function ConfigMethods:CREATE()
   local part = self.grid.partition[1];
   local filename = "config."..tostring(part)..".in"
   local fh = io.open(filename,"w");
   fh:write("\n! ------ BEGIN [",filename,"] file generated by <luacfg>\n\n");
   writeGRID(fh,self.grid);
   writeFDTD(fh,self.fdtd);
   writeBOUND(fh,self.bound);
   for _, src in ipairs(self.src) do writeSRC(fh,src) end
   for _, mat in ipairs(self.mat) do writeMAT(fh,mat) end
   for _, diag in ipairs(self.diag) do writeDIAG(fh,diag) end
   fh:write("! ------ END [",filename,"] \n\n");
   fh:close();
end


---------------------------------------------------------------------------
-- FINISH
---------------------------------------------------------------------------

-- move things to global namespace

for k,v in pairs(_M) do _G[k] = v end
for k,v in pairs(geo) do _G[k] = v end