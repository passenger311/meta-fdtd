
cfg = CONFIG{scenes=true}  -- scenes takes either true or false 

funcs = require "funcs"
round = funcs.round

dofile("scale.lua")

--- set time steps 

ncyc        = 20000


--- setup source / diagnostic planes for ffts

pulse = { 
   shape   = "Sech",
   width   = 200,
   offset  = 0,                    
   attack  = 10000,        
   sustain = 0, 
   decay   = 10000,   
}

planewave = { phi=00, theta=90, psi=00, nrefr=1 }


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
      amplitude = 5760000,
      pulse = pulse,
      planewave = planewave
   },
   REG{
      POINT{ 
	 { 10, ":", 0, 1 } 
      } 
   },
}


--- Definition of a three level system with E1 <= E2 <= E3
--- dipole matrix elements in units of dx and rates in units of 1/dt

cfg:MAT{
   THREELVL{ 
      invlambda = { 0, invwavelength, invwavelength }, --- { f12, f13, f23 }; f12+f23=f13!
      gamma = { 1e-7, 1e-7, 1e-7 }, --- dephasing constants (broadening) { gamma12, gamma13, gamma23 } 
      sigma = { 0, 1e-8, 1.e-8 },   --- relaxation constants {lvl2 -> lvl1 , lvl3 -> lvl1 , lvl3 -> lvl2 }
      mx = { { 0,0 }, { 0,0 }, { 0,0 } },  
      --- dipole matrix in x-direction { {Re(mu12),Im(mu12)}, {Re(mu13),Im(mu13)}, {Re(mu23),Im(mu23)} }
      my = { { 0,0 }, { 0.01,0 }, { 0.01,0 } }, --- y-direction
      mz = { { 0,0 }, { 0,0 }, { 0,0 } }, --- z-direction
      densities = { 1, 0, 0 }, --- initial occupation of three lvl system {n1,n2,n3}
      n = 10 --- number of three level systems per grid cell
   },
   REG{
	POINT{
	  { 500, ":", 1, 1, 1 }
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
