cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

size_tfsf=4

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
   dx = { conv*1e-9/1e40, 1, 1, 1 }

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
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, 0, 0, 1, 0, 0, 1 }
         } 
      }
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
      temp = temp,
      kmax = kmax,
      numk = numk,
   },
   REG{
      POINT{
         { 0, 0, 0, ":", 1, 1, 1 }
      }
   },
--[[
   OUT{
      file = { "GPL", "dens" },
      type = { "N", "S", ".F." },
      time = { 0, ncycles, 100 },
      REG{
         POINT{
            { 0, 0, 1 }
         }
      },
   },
--]]
   on = false
}

-- Spectrum diagnosis module
cfg:DIAG{
   PSPEC{
      file = "fft_ref1",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { 0, 0, 0 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_ref2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { 0, 0, 0 }
      }
   },
   on = diagpspec
}

--- CREATE: config.<part>.in
cfg:CREATE()

