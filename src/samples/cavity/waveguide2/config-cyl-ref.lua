
cfg = CONFIG{scenes=true}

dofile("scale.lua")

imin = -rwaveguide-size_pad
imax = rwaveguide+size_pad
jmin = -rwaveguide-size_pad
jmax = rwaveguide+size_pad
kmin = 0
kmax = 2*hhwaveguide

print("Computational window without PML:  ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                        ", ncycles)
max_tstep = math.floor(n_max*math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:  ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 3,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { jmin-size_pml, jmax+size_pml },   -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml }    -- -"- in z direction

}


--- CREATE SCENE

scene1 = Scene{
   value=1. -- constant background permittivity
}
cylinder1 = Cylinder{
   at={0,0,hhwaveguide},
   radius = rwaveguide,
   height= 2*hhwaveguide+2*size_pml
}
scene1:add{
   cylinder1,
--   depth=1,
   value = n_waveguide^2
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid1 = Grid{
   from={-rwaveguide-2,-rwaveguide-2,-size_pml},                    -- lower rear left corner of scene-object
   to={rwaveguide+2,rwaveguide+2,2*hhwaveguide+size_pml}  -- upper front right corner of scene-object
   -- if no further input is given 'from' and 'to' correspond directly to grid points
}
-- coarse grid for .VTK-preview
grid2 = Grid{
   yee=false,                                         -- one permittivity component only
   from={-rwaveguide-2,-rwaveguide-2,0},                    -- lower rear left corner of scene-object
   to={rwaveguide+2,rwaveguide+2,2*hhwaveguide}, -- upper front right corner of scene-object
   cells={50,50,200}                                   -- number of cells in each direction to grid scene
}
-- waveguide grid in injection plane
grid3 = Grid{ 
   from={-rwaveguide-1.5*resolution,-rwaveguide-1.5*resolution-1,kinj},
   to={rwaveguide+1.5*resolution,rwaveguide+1.5*resolution+1,kinj}
}


cfg:CREATE_GEO{      -- add the scene to the geometry
   "wg_ref",         -- filename: geo_"***".in
   scene=scene1,     -- scene to be added
   grid=grid1,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "wg_ref",         -- filename: preview_"***".vtk
   scene=scene1,     -- scene to be added
   grid=grid2        -- grid to be used
}

cfg:CREATE_GEO{      -- create the waveguide injection plane
   field_out,
   scene=scene1,
   grid=grid3
}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 LOAD_GEO{ "wg_ref" }
      },
      on = true
   },


   OUT{
      file = { "GPL", "ref_point_en_fft0" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, jmax, 3, kfft0, kfft0, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "ref_point_en_fft1" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft1, kfft1, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "ref_point_en_fft2" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft2, kfft2, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "ref_point_en_fft3" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft3, kfft3, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "ref_point_en_fft4" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft4, kfft4, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "ref_point_en_fft5" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{ 
            { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft5, kfft5, 1 }
         }
      }
   },

   OUT{
      file = { "VTK", "ref_xz_e" },
      type = { "E", "N" },
      time = { 0, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "ref_yz_e" },
      type = { "E", "N" },
      time = { 0, ncycles, 200 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "ref_xy_back_e" },
      type = { "E", "N" },
      time = { 0, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, kfft0, kfft0, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "ref_xy_front_e" },
      type = { "E", "N" },
      time = { 0, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, kfft1, kfft1, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "ref_xy_front_far_e" },
      type = { "E", "N" },
      time = { 0, ncycles, 200 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, kfft2, kfft2, 1 }
         }
      },
      on = true
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
      amplitude = 1.0,
      pulse = { 
	 shape="Gaussian", 
	 width=math.floor(widthl*resolution*n_max+.5),
	 offset=math.floor(offsetl*resolution*n_max+.5),
	 attack=math.floor(attackl*resolution*n_max+.5),
	 sustain=math.floor(sustainl*resolution*n_max+.5),
	 decay=math.floor(decayl*resolution*n_max+.5)
      },
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=nrefr }
   },
   REG{
      LOAD{ field_inj }
   },
   on = true
}


--- PSPEC
cfg:DIAG{
   PSPEC{
      file = "fft_ref-0",
      time = { 1, ncycles, 1 },
--      reffile = "fft_ref-1",
      mode = "Eap",
      phasewrap = { 1, 0 },
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft0, kfft0, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-1",
      time = { 1, ncycles, 1 },
      mode = "Eap",
      phasewrap = { 1, 0 },
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3, kfft1, kfft1, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-2",
      time = { 1, ncycles, 1 },
      phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, 3,  kfft2, kfft2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-3",
      time = { 1, ncycles, 1 },
      phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide, kfft3, kfft3, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-4",
      time = { 1, ncycles, 1 },
--      reffile = "fft_ref-3",
       phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide,  kfft4, kfft4, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-5",
      time = { 1, ncycles, 1 },
--      reffile = "fft_ref-3",
       phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
         { -rwaveguide, rwaveguide, 3, -rwaveguide, rwaveguide,  kfft5, kfft5, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()

--- Write file for Matlab waveguide field evaluation

foutput = io.open("lua.matlab","w+")
foutput:write(n_max, "\n", inv_wavelength, "\n")
foutput:write(-rwaveguide-1.5*resolution, " ", rwaveguide+1.5*resolution, " ", 1, " ",
              -rwaveguide-1.5*resolution, " ", rwaveguide+1.5*resolution, " ", 1, " ",
              kinj, " ", kinj, " ", 1, "\n")
foutput:write("geo_",field_out,".in\n")
foutput:write(field_inj, "\n")
foutput:close()
