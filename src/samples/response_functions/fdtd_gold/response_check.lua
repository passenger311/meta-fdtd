cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = tfsf_inj-size_pad
imax = 100*(-tfsf_inj)+size_pad

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
   krange = { 0, 0 }    -- -"- in z direction

}


--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, ":", eps_bg, eps_bg, eps_bg }
         },
         BOX{
            { 0, imax+size_pml, 1, ":", eps_infDL, eps_infDL, eps_infDL  }
         }
      },
      on = true
   },

   OUT{
      file = { "GPL", "ey_1" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 50 },
      REG{
         BOX{
            { 2, 2, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ey_2" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 50 },
      REG{
         BOX{
            { 4, 4, 1 }
         }
      }
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
      planewave = { phi=0., theta=90.0, psi=0.0, nrefr=nrefr },
      config = { 1,0,0,0,0,0 }
   },
   REG{
      BOX{
         { tfsf_inj, -tfsf_inj, 1, tfsf_inj, -tfsf_inj, 1, tfsf_inj, -tfsf_inj, 1 }
      }
   },
   on = true
}

--- MAT Definition(s)

if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{
         BOX{
            { 0, imax+size_pml, 1  }
         }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{
         BOX{
            { 0, imax+size_pml, 1  }
         }
      }
   }
end
if (mat == 'silver' or mat == 'carbon') then
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL2/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL2/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl2            -- delta epsilon
      },
      REG{
         BOX{
            { 0, imax+size_pml, 1  }
         }
      }
   }
end


cfg:DIAG{
   PSPEC{
      file = "fft_1",
      time = { 0, ncycles, 64 },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { 2, 2, 1 }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft_2",
      reffile = "fft_1",
      time = { 0, ncycles, 64 },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { 4, 4, 1 }
      }
   }
}
--- CREATE: config.<part>.in
cfg:CREATE()
