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

scene_cap_inf = Scene{
   value = n_bg^2 -- constant background permittivity
}
scene_cap = Scene{
   value = 0 -- constant background permittivity
}
strip = Box{
   from = {-hwidth,-hheight,-2},
   to = {hwidth,hheight,2}
}

-- round edges by removing box edges and adding cylinders
strip1 = Box{
   from = {-hwidth,-hheight,-2},
   to = {-hwidth+redges,-hheight+redges,2}
}
strip = BinaryAndNot{strip,strip1}
strip1 = Box{
   from = {hwidth-redges,-hheight,-2},
   to = {hwidth,-hheight+redges,2}
}
strip = BinaryAndNot{strip,strip1}
strip1 = Box{
   from = {-hwidth,hheight-redges,-2},
   to = {-hwidth+redges,hheight,2}
}
strip = BinaryAndNot{strip,strip1}
strip1 = Box{
   from = {hwidth-redges,hheight-redges,-2},
   to = {hwidth,hheight,2}
}
strip = BinaryAndNot{strip,strip1}

cyl = Cylinder{
   at = {-hwidth+redges,-hheight+redges,0},
   radius = redges,
   height = 4
}
strip = BinaryOr{strip,cyl}
cyl = Cylinder{
   at = {-hwidth+redges,hheight-redges,0},
   radius = redges,
   height = 4
}
strip = BinaryOr{strip,cyl}
cyl = Cylinder{
   at = {hwidth-redges,-hheight+redges,0},
   radius = redges,
   height = 4
}
strip = BinaryOr{strip,cyl}
cyl = Cylinder{
   at = {hwidth-redges,hheight-redges,0},
   radius = redges,
   height = 4
}
strip = BinaryOr{strip,cyl}
-- finished rounding edges

-- push strips to desired positions
strip1 = Transform{
   strip,
   move = {0,-hdist-hheight,0}
}
strip = Transform{
   strip,
   move = {0,hdist+hheight,0}
}

-- add strips to geometry
scene_cap:add{
   strip,
   value = 1
}
scene_cap_inf:add{
   strip,
   value = eps_infDL
}
scene_cap:add{
   strip1,
   value = 1
}
scene_cap_inf:add{
   strip1,
   value = eps_infDL
}

-- define grids
grid_strip1 = Grid{
   from = {-hwidth-3,-hdist-2*hheight-3,-1},
   to = {hwidth+3,-hdist+3,1}
}
grid_strip2 = Grid{
   from = {-hwidth-3,hdist-3,-1},
   to = {hwidth+3,hdist+2*hheight+3,1}
}
grid_prev_strip = Grid{
   yee = false,
   from = {-hwidth-3,-hdist-2*hheight-3,-1},
   to = {hwidth+3,hdist+2*hheight+3,1}
}


-- create geometries
cfg:CREATE_GEO{      -- add the scene to the geometry
   "strip1",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_strip1,       -- grid to be used
   method="default", 
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
-- on=false          -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "strip1_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_strip1,       -- grid to be used
   method="default", 
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "strip2",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_strip2,       -- grid to be used
   method="default",
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
-- on=false          -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "strip2_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_strip2,       -- grid to be used
   method="default",
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_PREVIEW{      -- add the scene to the geometry
   "resonator_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_prev_strip,       -- grid to be used
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_PREVIEW{      -- add the scene to the geometry
   "resonator",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_strip,       -- grid to be used
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}

--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, -1, 1, 1, ":", eps_bg, eps_bg, eps_bg }
         },
         LOAD_GEO{ "strip1_inf" },
         LOAD_GEO{ "strip2_inf" },
      },
      on = true
   },

   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { 0, math.floor(ncycles/10), 20 },
      REG{
         BOX{
            { -hdist_tfsf_i+2, hdist_tfsf_i-2, 10, -hdist_tfsf_j+2, hdist_tfsf_j-2, 10, 0, 0, 1 }
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
         LOAD_GEO{ "strip1" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "strip1" }
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
         LOAD_GEO{ "strip1" }
      }
   }
end

if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "strip2" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "strip2" }
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
         LOAD_GEO{ "strip2" }
      }
   }
end

-- diagnostics DFT/FFT
dofile("diagmode.lua")
-- dofile("diagpspec.lua")

--- CREATE: config.<part>.in
cfg:CREATE()
