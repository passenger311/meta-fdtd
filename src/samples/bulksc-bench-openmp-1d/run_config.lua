cfg = CONFIG{scenes=true}
taskid1 = select(1,...)
taskid2 = select(2,...)
dofile("scale.lua")

imin = -size_tfsf-size_pad
imax = size_tfsf+size_pad
jmin = 0
jmax = 0
kmin = 0
kmax = 0

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)

--- GRID Definition

cfg:GRID{

   dim = 1,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = {1,ncycles},                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { jmin, jmax },   -- -"- in y direction
   krange = { kmin, kmax },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

}

cfg:CHECKPOINT{

   load = false,
   save = false,
   detail = 3

}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, 0, 0, 1, 0, 0, 1, ":", n_bg^2, n_bg^2, n_bg^2 }
         }
      },
      on = true
   },

   OUT{  
      file = { "GPL", "en_sum" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, 0, 0, 1, 0, 0, 1 }
         } 
      },
      on = false
   },

   OUT{
      file = { "GPL", "Ey" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{ 
         POINT{
            { 0, 0, 0 }
         }
      },    
      on = false
   },
   OUT{
      file = { "GPL", "Ey2" },
      type = { "Ey", "N", ".T." },
      time = { 0, ncycles, 100 },
      REG{ 
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, 0, 0, 1, 0, 0, 1 }
         }
      },    
      on = false
   },

   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { 0, ncycles, 100 },
      REG{ 
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, 0, 0, 1, 0, 0, 1 }
         }
      },    
      on = false
   },

}

--- BOUND Definition

cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },

   PML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFBOX{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl*resolution/dt+.5),
         offset=math.floor(poffsetl*resolution/dt+.5),
         attack=math.floor(pattackl*resolution/dt+.5),
         sustain=math.floor(psustainl*resolution/dt+.5),
         decay=math.floor(pdecayl*resolution/dt+.5)
      },
      planewave = { phi=0, theta=90, psi=ppsi, nrefr=n_bg },
      config = {1,1,0,0,0,0}
   },
   REG{
      BOX{
            { -size_tfsf, size_tfsf, 1, 0, 0, 1, 0, 0, 1 }
      }
   },
   on = probe_on
}


cfg:MAT{
   BULKSC{
      gammap = gammap,
      M = M,
      egap = egap,
      me = me,
      mh = mh,
      N0 = N0,
      gammanr = gammanr,
      pump = 0,
      temp = temp,
      kmax = k_max,
      numk = taskid1,
   },
   REG{
      BOX{
         { -size_SC, size_SC, 1, 0, 0, 1, 0, 0, 1, ":", 1, 1, 1 }
      }
   },
   OUT{
      file = { "GPL", "carrier_dens" },
      type = { "N", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         POINT{
            { 0, 0, 0 }
         }
      },
    on = false
   },
   OUT{
      file = { "GPL", "carrier_pol" },
      type = { "P", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         POINT{
            { 0, 0, 0 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "carrier_dens" },
      type = { "N", "N" },
      time = { 0, ncycles, 1000 },
      REG{
         BOX{
            { -size_SC, size_SC, 1, 0, 0, 1, 0, 0, 1, ":", 1, 1, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "carrier_pol" },
      type = { "P", "N" },
      time = { 0, ncycles, 1000 },
      REG{
         BOX{
            { -size_SC, size_SC, 1, 0, 0, 1, 0, 0, 1, ":", 1, 1, 1 }
         }
      },
      on = false
   },
   on = true
}

-- Spectrum diagnosis module
--if (diagpspec) then
--  dofile("diagpspec.lua")
--end

--- CREATE: config.<part>.in
cfg:CREATE()

