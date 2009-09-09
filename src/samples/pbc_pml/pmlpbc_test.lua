cfg = CONFIG{scenes=true}
dofile("scale.lua")


imin = -hdist_ntff_ij
imax = hdist_ntff_ij
jmin = -hdist_ntff_k
jmax = hdist_ntff_k
kmin = 0
kmax = 0

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(n_max*math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 2,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { jmin, jmax },   -- -"- in y direction
   krange = { kmin, kmax }    -- -"- in z direction

}

--- FDTD Definition
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-size_pml-1, kmax+size_pml+1, 1, ":", n_bg^2, n_bg^2, n_bg^2 }
         },
         BOX{
            {-2,2,1,-2,2,1,0,0,1, ':', eps_infDL,eps_infDL,eps_infDL}
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "ex_liney" },
      type = { "Ex", "N" },
      time = { 0, math.floor(ncycles), 100 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, 0, 0, 1 }
         }
      }
   },
   OUT{
      file = { "VTK", "cap_xy_e" },
      type = { "E", "N" },
      time = { 100, ncycles, 100 },
      REG{
         BOX{
            { imin-size_pml, imax+size_pml, 1, jmin, jmax, 1, 0, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "en_point" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{
            { imax-1, imax-1, 1, jmax-1, jmax-1, 1, 0, 0, 1 }
         }
      }
   },


}

--- BOUND Definition

cfg:BOUND{

   config = { 1, 1, 3, 3, 1, 1 },

   PML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444,
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
         decay=math.floor(decayl*resolution*n_max+.5),
      },
      config = {0,0,1,0,0,0},
      planewave = { phi=90, theta=90.0, psi=0.0, nrefr=nrefr }
   },
   REG{
      BOX{
         {-hdist_tfsf_ij,hdist_tfsf_ij,1,-hdist_tfsf_k,hdist_tfsf_k,1,0,0,1}
      }
   },
   on = true
}

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
         decay=math.floor(decayl*resolution*n_max+.5),
         config = {0,0,0,1,0,0}
      },
      planewave = { phi=-90, theta=90.0, psi=180.0, nrefr=nrefr }
   },
   REG{
      BOX{
         {-hdist_tfsf_ij,hdist_tfsf_ij,1,-hdist_tfsf_k,hdist_tfsf_k,1,0,0,1}
      }
   },
   on = false
}


--- MAT Definition

if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         BOX{
            {-2,2,1,-2,2,1,0,0,1, ':', 1,1,1}
         }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         BOX{
            {-2,2,1,-2,2,1,0,0,1, ':', 1,1,1}
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
      REG{                                -- region where material is defined
         BOX{
            {-2,2,1,-2,2,1,0,0,1, ':', 1,1,1}
         }
      }
   }
end


cfg:DIAG{
   PSPEC{
      file = "fft_ref-zabs",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_ij+2, hdist_tfsf_ij-2, 1, hdist_tfsf_k-2, hdist_tfsf_k-2, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
