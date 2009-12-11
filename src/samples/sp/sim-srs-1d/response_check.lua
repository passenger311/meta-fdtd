cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = 2*tfsf_inj-size_pad
imax = fft+size_pad

print("Computational window without PML:         ", imin, imax, 0, 0, 0, 0)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(n_max*(imax-imin)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 1,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { 0, 0 },   -- -"- in y direction
   krange = { 0, 0 },   -- -"- in z direction
   dx = { conv*1e-6, 1, 1, 1 }

}


--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, -1, 1, 1, -1 ,1 , 1, ":", eps_bg, eps_bg, eps_bg }
         },
      },
      on = true
   },

   OUT{
      file = { "GPL", "ey_fft_r" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 40 },
      REG{
         BOX{
            { 2*tfsf_inj, 2*tfsf_inj, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ey_fft_t" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 40 },
      REG{
         BOX{
            { fft, fft, 1 }
         }
      }
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
   TFSFINJ{
      invlambda = inv_wavelength,
      amplitude = ampl_s / math.sqrt( real_c_si^2 / ((1e-6 * real_dx)^3 * real_eps0_si ) * real_chi3R ),
      pulse = { 
         shape="Gaussian",
         width=math.floor(widthl*resolution*n_max+.5),
         offset=math.floor(offsetl*resolution*n_max+.5),
         attack=math.floor(attackl*resolution*n_max+.5),
         sustain=math.floor(sustainl*resolution*n_max+.5),
         decay=math.floor(decayl*resolution*n_max+.5)
      },
      planewave = { phi=0., theta=90.0, psi=0.0, nrefr=nrefr },
   },
   REG{
      BOX{
         { tfsf_inj, tfsf_inj, 1 }
      }
   },
   on = true
}


cfg:SRC{
   TFSFINJ{
      invlambda = inv_wavelength + real_dx * real_invlambdaR,
      amplitude = ampl_p / math.sqrt( real_c_si^2 / ( (1e-6 * real_dx)^3 * real_eps0_si ) * real_chi3R ),
      pulse = {
         shape="Gaussian",
         width=math.floor(widthl*resolution*n_max+.5),
         offset=math.floor(offsetl*resolution*n_max+.5),
         attack=math.floor(attackl*resolution*n_max+.5),
         sustain=math.floor(sustainl*resolution*n_max+.5),
         decay=math.floor(decayl*resolution*n_max+.5)
      },
      planewave = { phi=0., theta=90.0, psi=0.0, nrefr=nrefr },
   },
   REG{
      BOX{
         { tfsf_inj, tfsf_inj, 1 }
      }
   },
   on = true
}

--- MAT Definition(s)


cfg:MAT{
    CHI3{             
    invlambdar = real_dx * real_invlambdaR,
    gammar     = real_gammaR * real_dt,
    chi3r      = real_c_si^2 / ( (1e-6 * real_dx)^3 * real_eps0_si ) * real_chi3R,
    chi3k      = 0
   },
   REG{
      BOX{
         { 0, fft-2, 1 }
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "fft_t",
--      reffile = "fft_ref",
      time = { 0, ncycles, 8 },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { fft, fft, 1 } 
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft_r",
--      reffile = "fft_ref",
      time = { 0, ncycles, 8 },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { 2*tfsf_inj, 2*tfsf_inj, 1 }
      }
   }
}
--- CREATE: config.<part>.in
cfg:CREATE()
