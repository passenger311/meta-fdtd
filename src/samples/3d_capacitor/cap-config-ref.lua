cfg = CONFIG{scenes=true}
dofile("scale.lua")

hdist_ntff = 10
hdist_tfsf = 8

imin = -hdist_ntff-size_pad
imax = hdist_ntff+size_pad
jmin = -hdist_ntff-size_pad
jmax = hdist_ntff+size_pad
kmin = -hdist_ntff-size_pad
kmax = hdist_ntff+size_pad

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(n_max*math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 3,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { jmin-size_pml, jmax+size_pml },   -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml, 0 },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

}

--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-size_pml-1, kmax+size_pml+1, 1, ":", eps_bg, eps_bg, eps_bg }
         }
      },
      on = true
   }

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
      planewave = { phi=0, theta=0.0, psi=0.0, nrefr=nrefr }
   },
   REG{
      BOX{
         {-hdist_tfsf,hdist_tfsf,1,-hdist_tfsf,hdist_tfsf,1,-hdist_tfsf,hdist_tfsf,1}
      }
   },
   on = true
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-zabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=0 }
   },
   REG{
      BOX{
         { -hdist_tfsf+2, hdist_tfsf-2, 1, -hdist_tfsf+2, hdist_tfsf-2, 1, -hdist_tfsf+2, -hdist_tfsf+2, 1 }
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "fft_ref-Sinj",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=0 }   
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -hdist_tfsf+1, -hdist_tfsf+1, 1 }
      }
   }
}


cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft-ref",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { 0, 1, 1, 0, 1, 1, -hdist_tfsf+1, -hdist_tfsf+1, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "F"
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -hdist_tfsf+1, -hdist_tfsf+1, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
