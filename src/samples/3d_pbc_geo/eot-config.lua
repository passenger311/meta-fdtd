cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = -hdist_ntff_i
imax = hdist_ntff_i
jmin = -hdist_ntff_j
jmax = hdist_ntff_j
kmin = -hdist_ntff_k-size_pad
kmax = hdist_ntff_k+size_pad

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
   irange = { imin, imax },   -- range of computational window in x direction
   jrange = { jmin, jmax },   -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml }    -- -"- in z direction

}


--- CREATE SCENE

scene_cap_inf = Scene{
   value = n_bg^2 -- constant background permittivity
}
scene_cap = Scene{
   value = 0 -- constant background permittivity
}

cap_box1 = Box{
   from={-hperiod-4,-2*hperiod-4,-hsio2height},
   to={hperiod+4,2*hperiod+4,hsio2height}
}

cap_box2 = Box{
   from={-hperiod-4,-2*hperiod-4,hsio2height},
   to={hperiod+4,2*hperiod+4,hsio2height+2*hmetalheight}
}

cap_box3 = Box{
   from={-hperiod-4,-2*hperiod-4,-hsio2height-2*hholeheight},
   to={hperiod+4,2*hperiod+4,-hsio2height}
}
cap_box4 = Box{
   from={-hxhole,-hperiod-hyhole,-hsio2height-2*hholeheight-2},
   to={hxhole,-hperiod+hyhole,-hsio2height}
}
cap_box5 = Box{
   from={-hyhole,hperiod-hxhole,-hsio2height-2*hholeheight-2},
   to={hyhole,hperiod+hxhole,-hsio2height}
}
cap_box3 = BinaryAndNot{cap_box3,cap_box4}
cap_box3 = BinaryAndNot{cap_box3,cap_box5}

scene_cap_inf:add{
   cap_box1,
   value = eps_sio2
}
scene_cap_inf:add{
   cap_box2,
   value = eps_infDL
}
scene_cap_inf:add{
   cap_box3,
   value = eps_infDL
}
scene_cap:add{
   cap_box2,
   value = 1.
}
scene_cap:add{
   cap_box3,
   value = 1.
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_capmetal = Grid{
   from = {imin-2, jmin-2,hsio2height-2},                    -- lower rear left corner of scene-object
   to = {imax+2, jmax+2, hsio2height+2*hmetalheight+2},  -- upper front right corner of scene-object
}
grid_caphole = Grid{
   from = {imin-2, jmin-2,-hsio2height-2*hholeheight-2},                    -- lower rear left corner of scene-object
   to = {imax+2, jmax+2, -hsio2height+2},  -- upper front right corner of scene-object
}

-- coarse grid for .VTK-preview
grid_prev_capmetal = Grid{
   yee=false,
   from = {imin-2, jmin-2,hsio2height-2},                    -- lower rear left corner of scene-object
   to = {imax+2, jmax+2, hsio2height+2*hmetalheight+2},  -- upper front right corner of scene-object
}
grid_prev_caphole = Grid{
   yee=false,
   from = {imin-2, jmin-2,-hsio2height-2*hholeheight-2},                    -- lower rear left corner of scene-object
   to = {imax+2, jmax+2, -hsio2height+2},  -- upper front right corner of scene-object
}


cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap1",         -- filename: geo_"***".in
   scene=scene_cap,  -- scene to be added
   grid=grid_caphole,    -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false        -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap1_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_caphole,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false         -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap2",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_capmetal,   -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false        -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap2_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_capmetal,      -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false         -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap1",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_caphole,     -- grid to be used
--   on=false
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap2",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_capmetal,    -- grid to be used
--   on=false
}

--- FDTD Definition
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-1, imax+1, 1, jmin-1, jmax+1, 1, kmin-size_pml-1, kmax+size_pml+1, 1, ":", n_bg^2, n_bg^2, n_bg^2 }
         },
         BOX{
            { imin-1, imax+1, 1, jmin-1, jmax+1, 1, -hsio2height, hsio2height, 1, ":", eps_sio2, eps_sio2, eps_sio2 }
         },
         LOAD_GEO{ "cap1_inf" },
         LOAD_GEO{ "cap2_inf" },
      },
      on = true
   },

   OUT{
      file = { "GPL", "en_point" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{
            { imax-1, imax-1, 1, jmax-1, jmax-1, 1, kmax-1, kmax-1, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sz_point-zin" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{ 
            { -hdist_tfsf_i, hdist_tfsf_i, 2, -hdist_tfsf_j, hdist_tfsf_j, 2, -hdist_tfsf_k+2, -hdist_tfsf_k+2, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sz_point-z" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{ 
            { -hdist_tfsf_i, hdist_tfsf_i, 2, -hdist_tfsf_j, hdist_tfsf_j, 2, -hdist_tfsf_k-2, -hdist_tfsf_k-2, 1 }
         }
      }
   },

   OUT{
      file = { "GPL", "sz_point+z" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         BOX{ 
            { -hdist_tfsf_i, hdist_tfsf_i, 2, -hdist_tfsf_j, hdist_tfsf_j, 2, hdist_tfsf_k+2, hdist_tfsf_k+2, 1 }
         }
      }
   },

   OUT{
      file = { "VTK", "cap_xz_e" },
      type = { "E", "N" },
      time = { 200, 12000, 200 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin-size_pml, kmax+size_pml, 1 }
         }
      },
      on = true
   },

   OUT{
      file = { "VTK", "cap_yz_e" },
      type = { "E", "N" },
      time = { 200, 12000, 200 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin-size_pml, kmax+size_pml, 1 }
         }
      },
      on = true
   }
}

--- BOUND Definition

cfg:BOUND{

   config = { 3, 3, 3, 3, 1, 1 },

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
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=nrefr },
      config = {0,0,0,0,1,1}
   },
   REG{
      BOX{
         {-hdist_tfsf_i-1,hdist_tfsf_i+1,1,-hdist_tfsf_j-1,hdist_tfsf_j+1,1,-hdist_tfsf_k,hdist_tfsf_k,1}
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
         LOAD_GEO{ "cap1" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap1" }
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
         LOAD_GEO{ "cap1" }
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
         LOAD_GEO{ "cap2" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap2" }
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
         LOAD_GEO{ "cap2" }
      }
   }
end


--- include diagnostics

dofile("diag.lua")

--- CREATE: config.<part>.in
cfg:CREATE()
