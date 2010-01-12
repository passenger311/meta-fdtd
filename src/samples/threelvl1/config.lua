
cfg = CONFIG{scenes=true}  -- scenes takes either true or false 

funcs = require "funcs"
round = funcs.round

dofile("scale.lua")

--- set time steps 

ncyc        = 40000


--- setup source / diagnostic planes for ffts

pulse1 = { 
   shape   = "Sech",
   width   = tau,
   offset  = 0,                    
   attack  = 20000,        
   sustain = 0, 
   decay   = 20000,   
}

pulse2  = {
   shape = "Sech",
   width = tau,
   offset = 20000,
   attack = 20000,
   sustain = 0,
   decay = 20000,
}

pulse3 = {
  shape = "Sech",
  width = tau,
  offset = 40000,
  attack = 20000,
  sustain = 0,
  decay = 20000,
}

pulse4 = {
  shape = "Sech",
  width = tau,
  offset = 60000,
  attack = 20000,
  sustain = 0,
  decay = 20000,
}


planewave1 = { phi=00, theta=90, psi=00, nrefr=1 }
planewave2 = { phi=00, theta=90, psi=90, nrefr=1 }
planewave3 = { phi=180, theta=90, psi=00, nrefr=1 }
planewave4 = { phi=180, theta=90, psi=90, nrefr=1 }

--- GRID Definition

cpml = 11
imin = 0
imax = 200

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

iminb = imin0+20
eps1 = 1.
eps2 = 1.2

d1 = math.floor(1/invwavelength/eps1*0.25)
d2 = math.floor(1/invwavelength/eps2*0.25)

epsbox = {}
epsbox[1] = { iminb, iminb+d1-1, 1, ":", eps1, eps1, eps1 }
epsbox[2] = { iminb+d1, iminb+d1+d2-1, 1, ":", eps2, eps2, eps2 }
iminb = iminb+d1+d2

for i=1,24 do
table.insert(epsbox, {iminb,iminb+d1-1,1,":",eps1,eps1,eps1})
table.insert(epsbox, {iminb+d1,iminb+d1+d2-1,1,":",eps2,eps2,eps2}) 
iminb = iminb+d1+d2
end


cfg:FDTD{

   EPSILON{
     REG{
       BOX{
	{imin0-1,imax0+1,1,":",eps_bg,eps_bg,eps_bg}
      }
     } 
   },

   OUT{
      file = { "R", "efield" },
      type = { "E", "N" },
      time = { 0, ncyc, 100000 },
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
      amplitude = 0.5*ENHL,
      pulse = pulse1,
      planewave = planewave2,
      alpha = 0
   },
   REG{
      POINT{ 
	 { 10, ":", 0, 1 } 
      } 
   },
}

cfg:SRC{
   TFSFINJ{
      invlambda = invwavelength,
      amplitude = 0.5*ENHL,
      pulse = pulse1,
      planewave = planewave1,
      alpha = 90
   },
   REG{
      POINT{
	{ 10, ":", 0, 1 }
      }
   },
}

cfg:SRC{
   TFSFINJ{
      invlambda = invwavelength,
      amplitude = 0*ENHL,
      pulse = pulse3,
      planewave = planewave1,
      alpha = 0
   },
   REG{
      POINT{
         { 10, ":", 0, 1 }
      }
   },
}

cfg:SRC{
   TFSFINJ{
      invlambda = invwavelength,
      amplitude = 0*ENHL,
      pulse = pulse4,
      planewave = planewave1,
      alpha = 0
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
      invlambda = { 0*invwavelength, invwavelength, invwavelength }, --- { f12, f13, f23 }; f12+f23=f13!
      gamma = { 0, 0, 0 }, --- dephasing constants (broadening) { gamma12, gamma13, gamma23 } 
      sigma = { 0, 0, 0 },   --- relaxation constants {lvl2 -> lvl1 , lvl3 -> lvl1 , lvl3 -> lvl2 }
      mx = { { 0,0 }, { 0,0 }, { 0,0 } },  
      --- dipole matrix in x-direction { {Re(mu12),Im(mu12)}, {Re(mu13),Im(mu13)}, {Re(mu23),Im(mu23)} }
      my = { { 0,0 }, { mu,0 }, { 0,0 } }, --- y-direction
      mz = { { 0,0 }, { 0,-mu }, { 0,0 } }, --- z-direction
      densities = { 1, 0, 0 }, --- initial occupation of three lvl system {n1,n2,n3}
      LFE = 0,
      n = 1 --- number of three level systems per grid cell
   },
   REG{
	POINT{
	  { 100, ":" , 1, 1, 1 }
	}
	},
   OUT{
      file = { "R", "dens" },
      type = { "N", "N", ".F." },
      time = { 0, ncyc, 2 },
      REG{
	 POINT{ 
	    { 100 }   -- middle point of the structure  
	 }
      },
      on = true
   },
   on = true
   
}



--- CREATE: config.<part>.in
cfg:CREATE()
