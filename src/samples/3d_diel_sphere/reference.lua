
cfg = CONFIG{scenes=true}

dofile("scale.lua")


htfsf = 4

imin = -htfsf-size_pad
imax = htfsf+size_pad
jmin = -htfsf-size_pad
jmax = htfsf+size_pad
kmin = -htfsf-size_pad
kmax = htfsf+size_pad

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
   krange = { kmin-size_pml, kmax+size_pml }    -- -"- in z direction

}

--- CREATE SCENE

scene_sphere = Scene{
   value=n_bg^2 -- constant background permittivity
}

diec_sphere = Sphere{
   at={0,0,0},
   radius = rsphere
}

scene_sphere:add{
   diec_sphere,
   value = n_sphere^2
}


-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_sphere = Grid{
   from={-rsphere-3,-rsphere-3,-rsphere-3},  -- lower rear left corner of scene-object
   to={rsphere+3,rsphere+3,rsphere+3}  -- upper front right corner of scene-object
}

-- coarse grid for .VTK-preview
grid_prev_sphere = Grid{
   yee=false,                               -- one permittivity component only
   from={-rsphere-3,-rsphere-3,-rsphere-3}, -- lower rear left corner of scene-object
   to={rsphere+3,rsphere+3,rsphere+3}  -- upper front right corner of scene-object
}

cfg:CREATE_GEO{  
   "sphere",         -- filename: geo_"***".in
   scene=scene_sphere,     -- scene to be added
   grid=grid_sphere,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=false           -- will the scene be updated
}
cfg:CREATE_PREVIEW{
   "sphere",         -- filename: preview_"***".vtk
   scene=scene_sphere,     -- scene to be added
   grid=grid_prev_sphere,        -- grid to be used
   on=false
}

--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-size_pml-1, kmax+size_pml+1, 1, ":", eps_bg, n_bg^2, eps_bg }
         },
--      LOAD_GEO{ "sphere" }
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
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=nrefr }
   },
   REG{
      BOX{
         {-htfsf,htfsf,1,-htfsf,htfsf,1,-htfsf,htfsf,1}
      }
   },
   on = true
}


--- DIAG PSPEC

cfg:DIAG{
   PSPEC{
      file = "fft-ref",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -htfsf+2, -htfsf+2, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft1-ref",
      time = { 0, ncycles, 1 },
      mode = "EHT"
   },
   REG{
      BOX{
         { -1, 0, 1, -1, 0, 1, -htfsf+2, -htfsf+2, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft2-ref",
      time = { 0, ncycles, 1 },
      mode = "HT"
   },
   REG{
      BOX{
         { -1, 0, 1, -1, 0, 1, -htfsf+1, -htfsf+1, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
