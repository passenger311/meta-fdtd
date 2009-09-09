cfg = CONFIG{scenes=true}
dofile("scale.lua")

jmax = offset1 + offsetN

--- GRID Definition

cfg:GRID{
   dim = dim,                       -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin, imax },   -- range of computational window in x direction
   jrange = { jmin-size_pml, jmax+size_pml },   -- -"- in y direction
   krange = { kmin, kmax }   -- -"- in z direction
}



--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
         { imin-1, imax+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-1, kmax+1, 1, ":", eps_bg, eps_bg, eps_bg }
         }
      },
      on = true
   } 
}

--- BOUND Definition

cfg:BOUND{

   config = { 3, 3, 1, 1, 3, 3 },

   PML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFINJ{
      invlambda = inv_wavelength,
      amplitude = 1.0,
      pulse = { 
         shape="Gaussian",
         width=math.floor(widthl*resolution*n_max+.5),
         offset=math.floor(offsetl*resolution*n_max+.5),
         attack=math.floor(attackl*resolution*n_max+.5),
         sustain=math.floor(sustainl*resolution*n_max+.5),
         decay=math.floor(decayl*resolution*n_max+.5)
      },
      planewave = { phi=90.0, theta=90.0, psi=psi, nrefr=nrefr }
   },
   REG{
      BOX{
         { imin, imax, 1, jinj, jinj, 1, kmin, kmax, 1}
      }
   },
   on = true
}

cfg:DIAG{
   PSPEC{
      file = "fft_cyl_ref",
      time = { 0, ncycles, (ncycles+1)/128 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=psi }
   },
   REG{
      BOX{
         { 0, 0, 1, jfft1, jfft1, 1, 0, 0, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
