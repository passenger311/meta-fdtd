cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = -hdist_tfsf_k-size_pad
imax = 100*hdist_tfsf_k+size_pad

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


--- CREATE SCENE

scene_np_inf = Scene{
   value =n_bg^2
}
scene_np = Scene{
   value =0.
}
box1 = Box{
   from={-rnp,-3,-3},
   to = {imax,3,3},
}
scene_np_inf:add{
   box1,
   value = eps_infDL
}
scene_np:add{
   box1,
   value = 1   -- filling factor for metal!!!!!!!!!!!!
}


-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)

grid_np = Grid{
   from={-rnp-2,-1,-1},
   to={imax+size_pml+3,1,1}
}
grid_prev_np = Grid{
   yee = false,
   from={-rnp-2,-1,-1},
   to={imax+size_pml+3,1,1}
}


cfg:CREATE_GEO{
   "np",
   scene=scene_np,
   grid=grid_np
}
cfg:CREATE_PREVIEW{
   "np",
   scene=scene_np,
   grid=grid_prev_np,
   on=false
}
cfg:CREATE_GEO{
   "np_inf",
   scene=scene_np_inf,
   grid=grid_np
}
cfg:CREATE_PREVIEW{
   "np_inf",
   scene=scene_np_inf,
   grid=grid_prev_np,
   on=false
}


--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, -3,3, 1, -3,3, 1, ":", eps_bg, eps_bg, eps_bg }
         },
         LOAD_GEO{ "np_inf" },
      },
      on = true
   },

   OUT{
      file = { "GPL", "ey_1" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { -2,-2,1 ,0,0,1,0,0,1}
         }
      }
   },
   OUT{
      file = { "GPL", "ey_2" },
      type = { "Ey", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { 2, 2, 1 }
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
         {-hdist_tfsf_k,hdist_tfsf_k+1,1,-hdist_tfsf_k,hdist_tfsf_k,1,-hdist_tfsf_k,hdist_tfsf_k,1}
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
      REG{                                -- region where material is defined
         LOAD_GEO{ "np" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "np" }
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
         LOAD_GEO{ "np" }
      }
   }
end


cfg:DIAG{
   PSPEC{
      file = "fft_1",
      time = { 0, ncycles, 4 },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -2, -2, 1 }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft_2",
      reffile = "fft_1",
      time = { 0, ncycles, 4 },
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
--- CREATE: config.<part>.in
cfg:CREATE()
