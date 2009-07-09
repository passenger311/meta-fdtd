
cfg = CONFIG{scenes=true}  -- scences takes either true or false 

TASKID = select(1,...)

funcs = require "funcs"
round = funcs.round

dofile("scale.lua")
dofile("geo.lua")


--- calculate time steps 

optical_length = ( grid.length_in + grid.length_out + math.abs(grid.dist_inout)  ) * nref.si

ncyc        = 8*round(optical_length/dt)       -- number of time steps


--- setup source / diagnostic planes for ffts

infty = 10000000000

pulse = { 
   shape   = "Gaussian",
   width   = 1000,
   offset  = 0,                    
   attack  = 2000,        
   sustain = infty, 
   decay   = 0,   
}


--- pml cells

cpml = 11     -- in real units it will be real_dx times cpml 


--- source / diagnostic planes for ffts

kl = 20

--- setup grids

grid0 = Grid{
   from={imin0-1,jmin0-1},
   to={imax0+1,jmax0+1} 
}

cfg:CREATE_PREVIEW{ "eps", scene=scene_eps, grid=grid0, silent=false, on=true }
cfg:CREATE_GEO{ "eps", scene=scene_eps, grid=grid0, comps=3, silent=false, on=true }

cfg:CREATE_PREVIEW{ "met", scene=scene_met, grid=grid0, silent=false, on=true }
cfg:CREATE_GEO{ "met", scene=scene_met, grid=grid0, comps=3, silent=false, on=true }


--- fire up matlab mode calculator!

--cmd = "./luacfg wgmode2j.lua geo_inj1.in "..tostring(invwavelength).." "..tostring(betaeff).." 1" 

--print(cmd)
--os.execute(cmd)


--- GRID Definition

cfg:GRID{

   dim = 2,  
   partition = { 0, 1 },
   ncyc = ncyc,
   dt = dt, 
   irange = { imin0,imax0 },
   jrange = { jmin0,jmax0 },
   krange = { 0,0 }
}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{ LOAD_GEO{ "eps" } },
      on = true   
   },

   OUT{
      file = { "GPL", "point_e_mid" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, 0, 0 }   -- middle point of the structure  
	 }
      }
   },

   OUT{
      file = { "VTK", "slice_xy" },
      type = { "E", "N" },
      time = { 0, ncyc, 100 },
      REG{
	 BOX{ 
	    { imin, imax, 1, jmin, jmax, 1 }   -- middle point of the structure  
	 }
      }
   },

}

--- BOUND Definition (pml boundary condition)

cfg:BOUND{
   config = { 1, 1, 1, 1, 1, 1 },   -- 1 pml 0 means not
   PML{
      cells = 11,                   -- no of pml cells
      pot = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s) (source)


xoff = 30
x0 = round(- grid.hwidth_cpl - grid.length_in + xoff)
y0 = round( grid.offset_in - grid.width_in/2)
y1 = y0+round(grid.width_in)


cfg:SRC{
   TFSFINJ{                             -- total field scattered field injection 
      invlambda = invwavelength,        -- inverse of wavelength (frequency)
      amplitude = 1.,                   -- amplitude of the source
      pulse = pulse,
      planewave = { phi=0, theta=0., psi=0, nrefr=nref.sio2 }
   },
   REG{
      BOX{ { x0, x0, 1, y0, y1, 1, ":", 1., -nref.sio2 } } 
   },
   on = true
}


--cfg:MAT{
--   PEC{},
--   DRUDE{
--      invlambdapl = invlambdapl,
--      gammapl = gammapl,
--      order = 2
--   },
--   REG{
--      LOAD_GEO{ "met" }
--   }
--}


--- CREATE: config.<part>.in
cfg:CREATE()