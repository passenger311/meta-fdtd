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
 

----------------- definitions


function Config()
   local tab = {}
   tab.MAT = {} -- subtables
   tab.SRC = {}
   tab.DIAG = {}
   return setmetatable(tab,ConfigMethods)
end


-- Grid Block definition

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

-- Fdtd Block definition

function ConfigMethods:FDTD(parms)
   self.FDTD = { block="FDTD" }
   for i,v in ipairs(parms) do self.FDTD[i] = parms[i] end
end

-- Epsilon Sub-Block definition

function EPSILON(parms)
   local EPSILON = { block="EPSILON" }
   for i,v in ipairs(parms) do EPSILON[i] = parms[i] end 
   return EPSILON
end

-- Bound Block definition

function ConfigMethods:BOUND(parms)
   self.BOUND = { block = "BOUND" }
   self.BOUND.config = parms.config or { 0,0,0,0,0,0 }
   for i = 1,6 do 
      if not self.BOUND.config[i] then self.BOUND.config[i] = 0 end 
   end
   for i,v in ipairs(parms) do self.BOUND[i] = parms[i] end
end

-- Pml Sub-Block definition

function PML(parms)
   local PML = { block = "PML" }
   PML.cells = parms.cells or 11
   PML.pot = parms.pot or 3.2
   PML.sigma = parms.sigma or 1.94444444444444
   PML.kappa = parms.kappa or 1.1
   for i,v in ipairs(parms) do PML[i] = parms[i] end
   return PML
end

-- Reg Sub-Block definition

function REG(parms) 
   local REG = { block = "REG" }
   REG.mask = parms.mask
   REG.list = parms.list
   REG.auto = parms.auto
   for i,v in ipairs(parms) do REG[i] = parms[i] end
   return REG
end

-- Load Sub-Block definition

function LOAD(parms) 
   local LOAD = { block = "LOAD" }
   for i,v in ipairs(parms) do LOAD[i] = parms[i] end
   return LOAD
end

-- Geo Sub-Block definition

function GEO(parms) 
   local GEO = { block = "GEO" }
   GEO.scene = parms.scene
   GEO.grid = parms.grid
   GEO.method = parms.method
   return GEO
end

-- Point Sub-Block definition

function POINT(parms) 
   local POINT = { block = "POINT" }
   for i,v in ipairs(parms) do POINT[i] = parms[i] end
   return POINT
end

-- Box Sub-Block definition

function BOX(parms) 
   local BOX = { block = "BOX" }
   for i,v in ipairs(parms) do BOX[i] = parms[i] end
   return BOX
end

-- Out Sub-Block definition

function OUT(parms) 
   local OUT = { block = "OUT" }
   OUT.file = parms.file or { "VTK", "?" }
   OUT.type = parms.type or { "?", "N", ".T." }
   OUT.type[2] = OUT.type[2] or "N"
   OUT.type[3] = OUT.type[3] or ".T."
   OUT.time = parms.time or { 0,-1, 1 }  
   for i,v in ipairs(parms) do OUT[i] = parms[i] end
   return OUT
end

-- Mat Block definition

function ConfigMethods:MAT(parms)
   local MAT = { block = "MAT" }
   for i,v in ipairs(parms) do MAT[i] = parms[i] end
   table.insert(self.MAT, MAT)   
end

-- Drude Sub-Block definition

function DRUDE(parms) 
   local DRUDE = { block = "DRUDE" }
   DRUDE.invlambdapl = parms.invlambdapl
   DRUDE.gammapl = parms.gammapl or 0
   DRUDE.order = parms.order or 1
   return DRUDE
end


-- Src Block definition

function ConfigMethods:SRC(parms)
   local SRC = { block = "SRC" }
   for i,v in ipairs(parms) do SRC[i] = parms[i] end
   table.insert(self.SRC, SRC)   
end




----------------- end of definitions; processing follows


local function writeGRID(fh,GRID)
   assert(GRID and GRID.block == "GRID", "Expected Grid{}")
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

local function writeBOX(fh,BOX)
   assert(BOX and BOX.block == "BOX", "Expected Box{}")
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

local function writePOINT(fh,POINT)
   assert(POINT and POINT.block == "POINT", "Expected Point{}")
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

local function writeLOAD(fh,LOAD)
   assert(LOAD and LOAD.block == "LOAD", "Expected Load{}")
   if LOAD.on == false then return end  
   fh:write("      (LOAD\n")
   fh:write("\t",LOAD.file,"\t! file to load\n")
   fh:write("      )LOAD\n")
end

local function writeREG(fh,REG)
   assert(REG and REG.block == "REG", "Expected Reg{}")
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

local function writeEPSILON(fh,EPSILON)
   assert(EPSILON and EPSILON.block == "EPSILON", "Expected Epsilon{}")
   if EPSILON.on == false then return end  
   fh:write("  (EPSILON\n")
   assert(EPSILON[1] and EPSILON[1].block == "REG","Epsilon{} must define Reg{}")
   writeREG(fh, EPSILON[1])
   fh:write("  )EPSILON\n")
end

local function writeOUT(fh,OUT)
   assert(OUT and OUT.block == "OUT", "Expected Out{}")
   assert(OUT[1] and OUT[1].block == "REG","Out{} must define Reg{}")
   if OUT.on == false then return end  
   fh:write("  (OUT\n")
   fh:write("    ",OUT.file[1]," ",OUT.file[2],"\t! file type and name\n"); 
   fh:write("    ",OUT.type[1]," ",OUT.type[2]," ",OUT.type[3],"\t! file type and name\n"); 
   fh:write("    ",OUT.time[1]," ",OUT.time[2]," ",OUT.time[3],"\t! time\n");
   writeREG(fh, OUT[1])
   fh:write("  )OUT\n")
end


local function writeFDTD(fh,FDTD)
   assert(FDTD and FDTD.block == "FDTD", "Expected Fdtd{}")
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

local function writePML(fh,PML)
   assert(PML and PML.block == "PML", "Expected Pml{}")
   fh:write("  (PML\n")
   fh:write("    ",PML.cells,"\t! # of cells\n");
   fh:write("    ",PML.pot,"\t! pot parameter\n");
   fh:write("    ",PML.sigma,"\t! sigma parameter\n");
   fh:write("    ",PML.kappa,"\t! kappa parameter\n");
   fh:write("  )PML\n")
end

local function writeBOUND(fh,BOUND)
   assert(BOUND and BOUND.block == "BOUND", "Expected Bound{}")
   fh:write("(BOUND\n")
   fh:write("  ",BOUND.config[1]," ",BOUND.config[2]," ",BOUND.config[3],
	    " ",BOUND.config[4]," ",BOUND.config[5]," ",BOUND.config[6],
	    "\t! boundary config\n"); 
   for i,v in ipairs(BOUND) do
      if v.block == "PML" then writePML(fh,v) end
   end
   fh:write(")BOUND\n\n")
end

local function writeSRC(fh,SRC)
   assert(SRC and SRC.block == "SRC", "Expected Src{}")
   if SRC.on == false then return end  
   local type = string.upper(SRC[1].block);
   fh:write("(SRC"..type.."\n")


   fh:write(")SRC"..type.."\n\n")
end

local writemat = {
   DRUDE = function(fh,DRUDE)
	      fh:write("  ", DRUDE.invlambdapl,"\t! invlambdapl [2 pi c]\n")
	      fh:write("  ", DRUDE.gammapl,"\t! gammapl (damping) [1/dt]\n")
	      fh:write("  ", DRUDE.order,"\t! order: 1 (J ode) or 2 (P ode)\n")
	   end
}

local function writeMAT(fh,MAT)
   assert(MAT and MAT.block == "MAT", "Expected Mat{}")
   if MAT.on == false then return end  
   assert( MAT[1] and MAT[2] and MAT[1].block and MAT[2].block, "Bad Mat{} structure") 
   local type = string.upper(MAT[1].block)
   fh:write("(MAT"..type.."\n")
   assert(writemat[type], "Mat{} type "..type.." does not exist")
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

local function writeDIAG(fh,DIAG)
   assert(DIAG and DIAG.block == "DIAG", "Expected Diag{}")
   if DIAG.on == false then return end  
   local type = string.upper(DIAG[1].block);
   fh:write("(DIAG"..type.."\n");
   

   fh:write(")DIAG"..type.."\n\n");
end



function ConfigMethods:CREATE()
   local part = self.GRID.partition[1];
   local fh = io.open("config."..tostring(part)..".in","w");
   fh:write("\n! ------ begin: file generated by <luacfg>\n\n");
   writeGRID(fh,self.GRID);
   writeFDTD(fh,self.FDTD);
   writeBOUND(fh,self.BOUND);
   for _, src in ipairs(self.SRC) do writeSRC(fh,src) end
   for _, mat in ipairs(self.MAT) do writeMAT(fh,mat) end
   for _, diag in ipairs(self.DIAG) do writeDIAG(fh,diag) end
   fh:write("! ------ end\n\n");
   fh:close();
end







----------------- move things into global namespace

for k,v in pairs(_M) do _G[k] = v end