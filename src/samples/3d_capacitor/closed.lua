cap_box3 = Box{
   at={-hwidthx+hheight,0,0},
   size={2*hheight,2*hwidthy,2*hdist+4*hheight}
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
   at={hwidthx-hheight,0,0},
   size={2*hheight,2*hwidthy,2*hdist+4*hheight}
}
scene_cap_inf:add{
   cap_box4,
   value = eps_infDL
}
scene_cap:add{
   cap_box4,
   value = 1.
}
cap_box5 = Box{
   at={0,-hwidthy+hheight,0,0},
   size={2*hwidthx,2*hheight,2*hdist+4*hheight}
}
scene_cap_inf:add{
   cap_box5,
   value = eps_infDL
}
scene_cap:add{
   cap_box5,
   value = 1.
}
cap_box6 = Box{
   at={0,hwidthy-hheight,0,0},
   size={2*hwidthx,2*hheight,2*hdist+4*hheight}
}
scene_cap_inf:add{
   cap_box6,
   value = eps_infDL
}
scene_cap:add{
   cap_box6,
   value = 1.
}

grid_cap3 = Grid{
   from = {-hwidthx+2*hheight+2,-hwidthy-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {-hwidthx-2,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}
grid_cap4 = Grid{
   from = {hwidthx-2*hheight-2,-hwidthy-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx+2,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}
grid_cap5 = Grid{
   from = {-hwidthx+2*hheight+3,-hwidthy+2*hheight+2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx-2*hheight-3,-hwidthy-2,hdist-3},  -- upper front right corner of scene-object
}
grid_cap6 = Grid{
   from = {-hwidthx+2*hheight+3,hwidthy-2*hheight-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx-2*hheight-3,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}

grid_prev_cap3 = Grid{
   yee = false,
   from = {-hwidthx+2*hheight+2,-hwidthy-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {-hwidthx-2,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}
grid_prev_cap4 = Grid{
   yee = false,
   from = {hwidthx-2*hheight-2,-hwidthy-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx+2,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}
grid_prev_cap5 = Grid{
   yee = false,
   from = {-hwidthx+2*hheight+3,-hwidthy+2*hheight+2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx-2*hheight-3,-hwidthy-2,hdist-3},  -- upper front right corner of scene-object
}
grid_prev_cap6 = Grid{
   yee = false,
   from = {-hwidthx+2*hheight+3,hwidthy-2*hheight-2,-hdist+3},                    -- lower rear left corner of scene-object
   to = {hwidthx-2*hheight-3,hwidthy+2,hdist-3},  -- upper front right corner of scene-object
}

cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap3",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap3,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap3_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap3,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap3",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap3,        -- grid to be used
   on=true
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap4",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap4,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap4_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap4,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap4",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap4,        -- grid to be used
   on=true
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap5",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap5,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap5_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap5,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap5",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap5,        -- grid to be used
   on=true
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap6",         -- filename: geo_"***".in
   scene=scene_cap,     -- scene to be added
   grid=grid_cap6,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_GEO{      -- add the scene to the geometry
   "cap6_inf",         -- filename: geo_"***".in
   scene=scene_cap_inf,     -- scene to be added
   grid=grid_cap6,       -- grid to be used
   method="default", -- ???
   comps=3,          -- number of components of permittivy (or permittivity + permeability)
   silent=false,     -- silent computation
   on=true           -- will the scene be updated
}
cfg:CREATE_PREVIEW{  -- create a preview of the scene
   "cap6",         -- filename: preview_"***".vtk
   scene=scene_cap,     -- scene to be added
   grid=grid_prev_cap6,        -- grid to be used
   on=true
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
if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap5" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap5" }
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
         LOAD_GEO{ "cap5" }
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
         LOAD_GEO{ "cap6" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "cap6" }
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
         LOAD_GEO{ "cap6" }
      }
   }
end

table.insert(geo,LOAD_GEO{ "cap3_inf" })
table.insert(geo,LOAD_GEO{ "cap4_inf" })
table.insert(geo,LOAD_GEO{ "cap5_inf" })
table.insert(geo,LOAD_GEO{ "cap6_inf" })
