cfg = CONFIG{scenes=true}
dofile("scale.lua")



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


--- CREATE SCENE

scene_bg = Scene{
   value=n_bg^2 -- constant background permittivity
}
scene_cyl = Scene{
   value =0.
}
for i,v in ipairs(rcyl) do
   cyl = Cylinder{
      at={icyl[i],jcyl[i],kcyl[i]},
      radius=rcyl[i],
      height=2*hhcyl[i]+1
   }
   scene_bg:add{
      cyl,
      value = eps_infDL
   }
   scene_cyl:add{
      cyl,
      value = 1   -- filling factor for metal
   }
end

table.sort(rcyl); table.sort(hhcyl); table.sort(icyl); table.sort(jcyl); table.sort(kcyl)
maxntab = table.maxn(rcyl);


-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_cyl = Grid{
   from={icyl[1]-rcyl[maxntab]-2,jcyl[1]-rcyl[maxntab]-2,kcyl[1]-hhcyl[maxntab]-2},
   to={icyl[maxntab]+rcyl[maxntab]+2,jcyl[maxntab]+rcyl[maxntab]+2,kcyl[maxntab]+hhcyl[maxntab]+2}
}
-- coarse grid for .VTK-preview
grid_prev_cyl = Grid{
   yee=false,
   from={icyl[1]-rcyl[maxntab]-2,jcyl[1]-rcyl[maxntab]-2,kcyl[1]-hhcyl[maxntab]-2},
   to={icyl[maxntab]+rcyl[maxntab]+2,jcyl[maxntab]+rcyl[maxntab]+2,kcyl[maxntab]+hhcyl[maxntab]+2},
   cells = 100
}

cfg:CREATE_GEO{      -- add the scene to the geometry
   "bg",         -- filename: geo_"***".in
   scene=scene_bg,     -- scene to be added
   grid=grid_cyl,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "bg",         -- filename: preview_"***".vtk
   scene=scene_bg,     -- scene to be added
   grid=grid_prev_cyl,        -- grid to be used
   on=true

}cfg:CREATE_GEO{
   "cyl",
   scene=scene_cyl,
   grid=grid_cyl
}
cfg:CREATE_PREVIEW{
   "cyl",
   scene=scene_cyl,
   grid=grid_prev_cyl,
   on=true
}


--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
         { imin-1, imax+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-1, kmax+1, 1, ":", eps_bg, eps_bg, eps_bg }
         },
         LOAD_GEO{ "bg" }
      },
      on = true
   },

   OUT{
      file = { "GPL", "en_point_cyl_inj" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 3, jinj, jinj, 1, kmin, kmax, dk }
         }
      }
   }, 

   OUT{
      file = { "GPL", "en_point_cyl_fft1" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 3, jfft1, jfft1, 1, kmin, kmax, dk }
         }
      }
   }, 
   
   OUT{
      file = { "GPL", "en_point_cyl_fft2" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 3, jfft2, jfft2, 1, kmin, kmax, dk }
         }
      }
   }, 

   OUT{
      file = { "VTK", "cyl_xz_inj_e" },
      type = { "E", "N" },
      time = { 200, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jinj, jinj, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },
   
   OUT{
      file = { "VTK", "cyl_xz_fft1_e" },
      type = { "E", "N" },
      time = { 200, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jfft1, jfft1, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },
   
   OUT{
      file = { "VTK", "cyl_xz_fft2_e" },
      type = { "E", "N" },
      time = { 200, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jfft2, jfft2, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },
      
   OUT{
      file = { "VTK", "cyl_yz_e" },
      type = { "E", "N" },
      time = { 200, ncycles, 200 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "cyl_xy_e" },
      type = { "E", "N" },
      time = { 200, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, 0, 1 }
         }
      },
      on = true
   },
   
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

--- MAT Definition(s)

cfg:MAT{                               -- define material with frequency/time-dependent answer
   DRUDE{                              -- define a Drude material
      invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
      gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
      order = 2
   },
   REG{                                -- region where material is defined
      LOAD_GEO{ "cyl" }
   }
}
cfg:MAT{
   LORENTZ{                            -- define a Lorentzian material
      invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
      gammal = conv*real_gammaL/frequ_factor,       -- resonance width
      deltaepsl = deltaepsl            -- delta epsilon
   },
   REG{                                -- region where material is defined
      LOAD_GEO{ "cyl" }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_cyl_R",
      reffile = "fft_cyl_ref",
      time = { 0, ncycles, (ncycles+1)/128 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=psi }
   },
   REG{
      BOX{
         { imin, imax, 3, jfft1, jfft1, 1, kmin, kmax, dk }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft_cyl_T",
      reffile = "fft_cyl_ref",
      time = { 0, ncycles, (ncycles+1)/128 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=psi }
   },
   REG{
      BOX{
         { imin, imax, 3, jfft2, jfft2, 1, kmin, kmax, dk }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
