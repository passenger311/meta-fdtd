cfg = CONFIG{scenes=true}
taskid = select(1,...)

dofile("scale.lua")

imin = -hdist_ntff_i-size_pad
imax = hdist_ntff_i+size_pad
jmin = -hdist_ntff_j-size_pad
jmax = hdist_ntff_j+size_pad
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
   jrange = { jmin-size_pml, jmax+size_pml },   -- -"- in y direction
   krange = { kmin, kmax },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

}


--- CREATE SCENE

scene_core_shell_inf = Scene{
   value = n_bg^2
}
scene_shell = Scene{
   value = 0
}
scene_core = Scene{
   value = 0
}
for i,v in ipairs(rnp) do
   sphere1 = Sphere{
      at={inp[i],jnp[i],knp[i]},
      radius=rnp[i]
   }
   sphere2 = Sphere{
      at={inp[i],jnp[i],knp[i]},
      radius=rnp[i]-shell[i]
   }
   shell = BinaryAndNot{sphere1,sphere2}
   scene_core_shell_inf:add{
      shell,
      value = eps_diel
   }
   scene_shell:add{
      shell,
      value = 1
   }
   scene_core_shell_inf:add{
      sphere2,
      value = eps_infDL
   }
   scene_core:add{
      sphere2,
      value = 1   -- filling factor for metal!!!!!!!!!!!!
   }
end


-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
table.sort(rnp); table.sort(inp); table.sort(jnp); table.sort(knp)
maxntab = table.maxn(rnp);

grid_np = Grid{
   from={inp[1]-rnp[maxntab]-2,jnp[1]-rnp[maxntab]-2,-1},
   to={inp[maxntab]+rnp[maxntab]+2,jnp[maxntab]+rnp[maxntab]+2,1}
}
grid_prev_np = Grid{
   yee = false,
   from={inp[1]-rnp[maxntab]-2,jnp[1]-rnp[maxntab]-2,-1},
   to={inp[maxntab]+rnp[maxntab]+2,jnp[maxntab]+rnp[maxntab]+2,1}
}

-- create geometry
cfg:CREATE_GEO{
   "core",
   scene=scene_core,
   grid=grid_np,
--   on=false
}
cfg:CREATE_PREVIEW{
   "core",
   scene=scene_core,
   grid=grid_prev_np,
   on=false
}

cfg:CREATE_GEO{
   "shell",
   scene=scene_shell,
   grid=grid_np,
--   on=false
}
cfg:CREATE_PREVIEW{
   "shell",
   scene=scene_shell,
   grid=grid_prev_np,
   on=false
}
cfg:CREATE_GEO{
   "core_shell_inf",
   scene=scene_core_shell_inf,
   grid=grid_np,
--   on=false
}
cfg:CREATE_PREVIEW{
   "core_shell_inf",
   scene=scene_core_shell_inf,
   grid=grid_prev_np,
--   on=false
}

--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, -1, 1, 1, ":", eps_bg, eps_bg, eps_bg }
         },
         LOAD_GEO{ "core_shell_inf" }
      },
      on = true
   },

   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { 0, math.floor(ncycles/10), 20 },
      REG{
         BOX{
            { -hdist_tfsf_i+2, hdist_tfsf_i-2, 4, -hdist_tfsf_j+2, hdist_tfsf_j-2, 4, 0, 0, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "en_point" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 5 },
      REG{
         BOX{
            { imax-1, imax-1, 1, jmax-1, jmax-1, 1, 0, 0, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sy_point-yin" },
      type = { "Sy", "S", ".F." },
      time = { 0, ncycles, 5 },
      REG{
         BOX{ 
            { -hdist_tfsf_i+2, hdist_tfsf_i-2, 2, -hdist_tfsf_j+2, -hdist_tfsf_j+2, 1, 0, 0, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sy_point-y" },
      type = { "Sy", "S", ".F." },
      time = { 0, ncycles, 5 },
      REG{
         BOX{ 
            { -hdist_tfsf_i-2, hdist_tfsf_i+2, 2, -hdist_tfsf_j-2, -hdist_tfsf_j-2, 1, 0, 0, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sy_point+y" },
      type = { "Sy", "S", ".F." },
      time = { 0, ncycles, 5 },
      REG{
         BOX{ 
            { -hdist_tfsf_i-2, hdist_tfsf_i+2, 2, hdist_tfsf_j-2, hdist_tfsf_j-2, 1, 0, 0, 1 }
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
      planewave = { phi=90.0+5*(taskid-1), theta=90.0, psi=0.0, nrefr=nrefr }
   },
   REG{
      BOX{
         {-hdist_tfsf_i,hdist_tfsf_i,1,-hdist_tfsf_j,hdist_tfsf_j,1,-hdist_tfsf_j,hdist_tfsf_j,1}
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
         LOAD_GEO{ "core" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "core" }
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
         LOAD_GEO{ "core" }
      }
   }
end


-- diagnostics DFT/FFT
dofile("diagmode.lua")
-- dofile("diagpspec.lua")

--- CREATE: config.<part>.in
cfg:CREATE()
