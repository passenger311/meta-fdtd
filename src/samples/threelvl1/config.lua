
cfg = CONFIG{scenes=true}  -- scences takes either true or false 

funcs = require "funcs"
round = funcs.round

dofile("scale.lua")

--- set time steps 

ncyc        = 20000


--- setup source / diagnostic planes for ffts

pulse = { 
   shape   = "Gaussian",
   width   = 200,
   offset  = 0,                    
   attack  = 10000,        
   sustain = 0, 
   decay   = 10000,   
}

planewave = { phi=0, theta=90, psi=00, nrefr=1 }


--- GRID Definition

cpml = 11
imin = 0
imax = 1000

imin0 = imin - cpml
imax0 = imax + cpml

cfg:GRID{

   dim = 1,  
   partition = { 0, 1 },
   ncyc = ncyc,
   dt = dt, 
   irange = { imin0,imax0 },
   jrange = { 0,0 },
   krange = { 0,0 },
   dx = { real_dx*1e-6, 1, 1, 1 }   
}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
     REG{
      BOX{
        { imin0-1, imax0+1, 1, ":", eps_bg, eps_bg, eps_bg }
      }
     }
   },

   OUT{
      file = { "GPL", "efield" },
      type = { "E", "N" },
      time = { 0, ncyc, 1000 },
      REG{
	 BOX{ 
	    { imin, imax, 1 }   -- middle point of the structure  
	 }
      }
   },

}

--- BOUND Definition (pml boundary condition)

cfg:BOUND{
   config = { 1, 1, 1, 1, 1, 1 },   -- 1 pml 0 means not
   PML{
      cells = cpml,                   -- no of pml cells
      pot = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s) (source)


cfg:SRC{
   TFSFINJ{ 
      invlambda = invwavelength,
      amplitude = 1.,
      pulse = pulse,
      planewave = planewave
   },
   REG{
      POINT{ 
	 { 10, ":", 0, 1 } 
      } 
   },
}

cfg:MAT{
   THREELVL{ 
      invlambda = { invwavelength, invwavelength, 0 },
      gamma = { 0, 0, 0 },
      sigma = { 0, 0, 0 },
      mx = { { 0,0 }, { 0,0 }, { 0,0 } },
      my = { { 0.1,0 }, { 0.1,0 }, { 0,0 } },
      mz = { { 0,0 }, { 0,0 }, { 0,0 } },
      densities = { 1., 0., 0. },
      n = 10 
   },
   REG{
	POINT{
	  { 500 , ":", 1, 1, 1 }	
	}
   },
   OUT{
      file = { "GPL", "dens" },
      type = { "N", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 500 }   -- middle point of the structure  
	 }
      },
      on = true
   },
   on = true
   
}



--- CREATE: config.<part>.in
cfg:CREATE()