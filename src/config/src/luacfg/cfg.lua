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


local _G,print,pairs,ipairs,type,assert,setmetatable,table,string,io,tostring
   = 
   _G,print,pairs,ipairs,type,assert,setmetatable,table,string,io,tostring
local geo = require "geo"

module("cfg")

-- Config

ConfigMethods = {}
ConfigMethods.__index = ConfigMethods
 

---------------------------------------------------------------------------
-- BLOCK DEFINITIONS
---------------------------------------------------------------------------

function CONFIG()
   local tab = {}
   tab.mat = {} -- subtables
   tab.src = {}
   tab.diag = {}
   return setmetatable(tab,ConfigMethods)
end


-- GRID Config-Block definition

function ConfigMethods:GRID(parms)
   self.GRID = { block="GRID" }
   self.GRID.irange = parms.irange or { 0, 0 };
   self.GRID.jrange = parms.jrange or { 0, 0 };
   self.GRID.krange = parms.krange or { 0, 0 };
   self.GRID.dim =  parms.dim or 1;
   self.GRID.partition = parms.partition or { 0,1 };
   self.GRID.ncyc = parms.ncyc or 100;
   self.GRID.dt = parms.ncyc or 0.9999;
end

-- FDTD Config-Block definition

function ConfigMethods:FDTD(parms)
   self.FDTD = { block="FDTD" }
   for i,v in ipairs(parms) do self.FDTD[i] = parms[i] end
end

-- EPSILON Sub-Block definition

function EPSILON(parms)
   local EPSILON = { block="EPSILON" }
   for i,v in ipairs(parms) do EPSILON[i] = parms[i] end 
   return EPSILON
end

-- BOUND Config-Block definition

function ConfigMethods:BOUND(parms)
   self.BOUND = { block = "BOUND" }
   self.BOUND.config = parms.config or { 0,0,0,0,0,0 }
   for i = 1,6 do 
      if not self.BOUND.config[i] then self.BOUND.config[i] = 0 end 
   end
   for i,v in ipairs(parms) do self.BOUND[i] = parms[i] end
end

-- PML Sub-Block definition

function PML(parms)
   local PML = { block = "PML" }
   PML.cells = parms.cells or 11
   PML.pot = parms.pot or 3.2
   PML.sigma = parms.sigma or 1.94444444444444
   PML.kappa = parms.kappa or 1.1
   for i,v in ipairs(parms) do PML[i] = parms[i] end
   return PML
end

-- REG Sub-Block definition

function REG(parms) 
   local REG = { block = "REG" }
   REG.mask = parms.mask
   REG.list = parms.list
   REG.auto = parms.auto
   for i,v in ipairs(parms) do REG[i] = parms[i] end
   return REG
end

-- LOAD Sub-Block definition

function LOAD(parms) 
   local LOAD = { block = "LOAD" }
   for k,v in pairs(parms) do LOAD[k] = v end
   return LOAD
end

-- (NEW!) GEO Sub-Block definition

function GEO(parms) 
   local GEO = { block = "GEO" }
   for k,v in pairs(parms) do GEO[k] = v end
   return GEO
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

-- DIAG Config-Block definition

function ConfigMethods:DIAG(parms)
   local DIAG = { block = "DIAG" }
   for k,v in pairs(parms) do DIAG[k] = v end
   table.insert(self.diag, DIAG)   
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
   fh:write("      (BOX\n")
   for b=1, #BOX, 2 do
      fh:write("\t")
      coord = BOX[b]
      fillf = BOX[b+1]
      for i=1,#coord do
	 fh:write(coord[i]," ");
      end
      if fillf and #fillf > 0 then
	 fh:write(": ");
	 for i=1,#fillf do
	    fh:write(fillf[i]," ");
	 end
      end
      fh:write("\n");
   end
   fh:write("      )BOX\n")
end

-- POINT Sub-Block write

local function writePOINT(fh,POINT)
   assert(POINT and POINT.block == "POINT", "Expected POINT{}")
   if POINT.on == false then return end  
   fh:write("      (POINT\n")
   for p=1, #POINT, 2 do
      fh:write("\t")
      coord = POINT[p]
      fillf = POINT[p+1]
      for i=1,#coord do
	 fh:write(coord[i]," ");
      end
      if fillf and #fillf > 0 then
	 fh:write(": ");
	 for i=1,#fillf do
	    fh:write(fillf[i]," ");
	 end
      end
      fh:write("\n");
   end
   fh:write("      )POINT\n")
end

-- LOAD Sub-Block write

local function writeLOAD(fh,LOAD)
   assert(LOAD and LOAD.block == "LOAD", "Expected LOAD{}")
   if LOAD.on == false then return end  
   fh:write("      (LOAD\n")
   fh:write("\t",LOAD.file,"\t! file to load\n")
   fh:write("      )LOAD\n")
end

-- GEO Sub-Block write

local geonum = 0

local function writeGEO(fh,GEO)
   assert(GEO and GEO.block == "GEO", "Expected GEO{}")
   if GEO.on == false then return end  
   geonum = geonum + 1
   GEO.file = "config_geo"..tostring(geonum)..".in";
   geo_fh = geo.FileIN(GEO.file);
   assert(GEO.scene and GEO.grid,"GEO{} must define <grid> and <scene>") 
   print("processing "..GEO.file.." ... ")
   GEO.grid:write{geo_fh,GEO.scene}
   print("DONE.\n")
   fh:write("      (LOAD\n")
   fh:write("\t",GEO.file,"\t! file to load\n")
   fh:write("      )LOAD\n")
end

-- REG Sub-Block write

local function writeREG(fh,REG)
   assert(REG and REG.block == "REG", "Expected REG{}")
   if REG.on == false then return end  
   fh:write("    (REG\n")
   for i,v in ipairs(REG) do
      if v.block == "BOX" then writeBOX(fh,v) end
      if v.block == "POINT" then writePOINT(fh,v) end
      if v.block == "LOAD" then writeLOAD(fh,v) end
      if v.block == "GEO" then writeGEO(fh,v) end
   end
   if REG.auto then fh:write("      AUTO\t! auto loop mode\n") end
   if REG.mask then fh:write("      MASK\t! mask loop mode\n") end
   if REG.list then fh:write("      LIST\t! list loop mode\n") end
   fh:write("    )REG\n")
end

-- EPSILON Sub-Block write

local function writeEPSILON(fh,EPSILON)
   assert(EPSILON and EPSILON.block == "EPSILON", "Expected EPSILON{}")
   if EPSILON.on == false then return end  
   fh:write("  (EPSILON\n")
   assert(EPSILON[1] and EPSILON[1].block == "REG","EPSILON{} must define REG{}")
   writeREG(fh, EPSILON[1])
   fh:write("  )EPSILON\n")
end

-- OUT Sub-Block write

local function writeOUT(fh,OUT)
   assert(OUT and OUT.block == "OUT", "Expected OUT{}")
   assert(OUT[1] and OUT[1].block == "REG","OUT{} must define REG{}")
   if OUT.on == false then return end  
   fh:write("  (OUT\n")
   fh:write("    ",OUT.file[1]," ",OUT.file[2],"\t! file type and name\n"); 
   fh:write("    ",OUT.type[1]," ",OUT.type[2]," ",OUT.type[3],"\t! file type and name\n"); 
   fh:write("    ",OUT.time[1]," ",OUT.time[2]," ",OUT.time[3],"\t! time\n");
   writeREG(fh, OUT[1])
   fh:write("  )OUT\n")
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
   fh:write("  (PML\n")
   fh:write("    ",PML.cells,"\t! # of cells\n");
   fh:write("    ",PML.pot,"\t! pot parameter\n");
   fh:write("    ",PML.sigma,"\t! sigma parameter\n");
   fh:write("    ",PML.kappa,"\t! kappa parameter\n");
   fh:write("  )PML\n")
end

-- BOUND Config-Block write

local function writeBOUND(fh,BOUND)
   assert(BOUND and BOUND.block == "BOUND", "Expected BOUND{}")
   fh:write("(BOUND\n")
   fh:write("  ",BOUND.config[1]," ",BOUND.config[2]," ",BOUND.config[3],
	    " ",BOUND.config[4]," ",BOUND.config[5]," ",BOUND.config[6],
	    "\t! boundary config\n"); 
   for i,v in ipairs(BOUND) do
      if v.block == "PML" then writePML(fh,v) end
   end
   fh:write(")BOUND\n\n")
end

-- (VARIOUS) SRC Sub-Block write

local writesrc = {

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

-- (VARIOUS) MAT Sub-Block write

local writemat = {
   DRUDE = function(fh,DRUDE)
	      fh:write("  ", DRUDE.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write("  ", DRUDE.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write("  ", DRUDE.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
   LHM = function(fh,LHM)
	      fh:write("  ", LHM.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write("  ", LHM.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write("  ", LHM.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end,
   LORENTZ = function(fh,LORENTZ)
	      fh:write("  ", LORENTZ.invlambdal,"\t! invlambdal [2 pi c]\n")
	      fh:write("  ", LORENTZ.gammal,"\t! gammal (damping) [1/dt]\n")
	      fh:write("  ", LORENTZ.deltaepsl,"\t! deltaepsl coupling coeff.\n")
	   end,
   DEBYE = function(fh,DEBYE)
	      fh:write("  ", DEBYE.taud,"\t! taud [dt]\n")
	      fh:write("  ", DEBYE.deltaepsd,"\t! deltaepsd coupling coeff.\n")
	   end,
   PEC = function(fh,PEC)
	 end,
   BLOCH = function(fh,BLOCH)
	      fh:write("  ", BLOCH.invlambdal,"\t! invlambdal [2 pi c]\n")
	      fh:write("  ", BLOCH.gammal,"\t! gammal (damping) [1/dt]\n")
	      fh:write("  ", 
		       "(",BLOCH.dipole[1][1],",",BLOCH.dipole[1][2],") ",
		       "(",BLOCH.dipole[2][1],",",BLOCH.dipole[2][2],") ",
		       "(",BLOCH.dipole[3][1],",",BLOCH.dipole[3][2],") ",
		       "\t! dipole matrix elem. []\n")
	      fh:write("  ", BLOCH.carrier[1]," ", BLOCH.carrier[2],
		       "\t! carrier numbers transp./initial  []\n")
	      fh:write("  ", BLOCH.gammanr,"\t! non-rad. recomb. [1/dt]\n")
	      fh:write("  ", BLOCH.pump,"\t! pump rate [1/dt]\n")
	      fh:write("  ", BLOCH.satmodel,"\t! sat.model 0=>(N-Ntr), 1=>Ntr*log(N/Ntr)\n")		       
	   end
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

-- (VARIOUS) DIAG Sub-Block write

local writediag = {

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
   local part = self.GRID.partition[1];
   local filename = "config."..tostring(part)..".in"
   local fh = io.open(filename,"w");
   fh:write("\n! ------ BEGIN [",filename,"] file generated by <luacfg>\n\n");
   writeGRID(fh,self.GRID);
   writeFDTD(fh,self.FDTD);
   writeBOUND(fh,self.BOUND);
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