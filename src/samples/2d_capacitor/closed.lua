cap_box3 = Box{
   from={-hwidth,-hdist-2*hheight,-3},
   to={-hwidth+2*hheight,hdist+2*hheight,3}
}
scene_cap_inf:add{
   cap_box3,
   value = eps_infDL
}
scene_cap:add{
   cap_box3,
   value = 1.
}

cap_box4 = Box{
   from={hwidth-2*hheight,-hdist-2*hheight,-3},
   to={hwidth,hdist+2*hheight,3}
}
scene_cap_inf:add{
   cap_box4,
   value = eps_infDL
}
scene_cap:add{
   cap_box4,
   value = 1.
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_cap3 = Grid{
   from = {-hwidth-2,-hdist+3,-2},                    -- lower rear left corner of scene-object
   to = {-hwidth+2*hheight+2,hdist-3,2},  -- upper front right corner of scene-object
}
grid_cap4 = Grid{
   from = {hwidth-2*hheight-2,-hdist+3,-2},                    -- lower rear left corner of scene-object
   to = {hwidth+2,hdist-3,2},  -- upper front right corner of scene-object
}
-- coarse grid for .VTK-preview
grid_prev_cap3 = Grid{
   yee = false,
   from = {-hwidth-2,-hdist+3,-2},                    -- lower rear left corner of scene-object
   to = {-hwidth+2*hheight+2,hdist-3,2},  -- upper front right corner of scene-object
}
grid_prev_cap4 = Grid{
   yee = false,
   from = {hwidth-2*hheight-2,-hdist+3,-2},                    -- lower rear left corner of scene-object
   to = {hwidth+2,hdist-3,2},  -- upper front right corner of scene-object
}

cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap3",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap3,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap3_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap3,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap3",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap3,        -- grid to be used
--   on=false
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap4",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap4,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap4_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap4,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
--   on=false          -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap4",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap4,        -- grid to be used
--   on=false
}

if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap3" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap3" }
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
         LOAD_GEO{ "cap3" }
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
         LOAD_GEO{ "cap4" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap4" }
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
         LOAD_GEO{ "cap4" }
      }
   }
end

table.insert(geo,LOAD_GEO{ "cap3_inf" })
table.insert(geo,LOAD_GEO{ "cap4_inf" })

