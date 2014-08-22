-- basic objects
box_metal_sub = Box{
  from = {imin-size_pml-1,jmin-size_pml-1,0},
  to = {imax+size_pml+size_pml+1,jmax+size_pml+1,-hmetal_sub-size_pml-1}
}
spacer = Box{
  from = {imin-size_pml-1,jmin-size_pml-1,hspacer},
  to = {imax+size_pml+1,jmax+size_pml+1,0}
}
-- box_metal_cyl = Box{
--   from = {-cyl_width,-cyl_length,-cyl_diameter-hspacer},
--   to = {cyl_width,cyl_length,-hspacer}
-- }

-- complex objects

metal_cyl = Cylinder{
  radius = cyl_radius,
  height = 2*cyl_length,
  at = {0, 0, 0}
}
metal_cyl = Transform{metal_cyl,angle=90,axis={1, 0, 0}}
scene_background = Scene{
   value = n_bg^2
}
metal_cyl = Transform{metal_cyl,move={0,0,cyl_radius+hspacer}}
metal_cyl = Transform{metal_cyl,move={-cyl_radius+math.floor(cyl_radius),0,0}}
scene_background:add{
   box_metal_sub,
   value = eps_infDL
}
scene_background:add{
   spacer,
   value = n_spacer^2
}
scene_background:add{
   metal_cyl,
   value = n_nanowire^2
}
scene_metal_sub = Scene{
   value = 0
}
scene_metal_sub:add{
   box_metal_sub,
   value = 1
}
--[[
scene_metal_cyl = Scene{
   value = 0
}
scene_metal_cyl:add{
   metal_cyl,
   value = 1
}
--]]
scene_4lvl = Scene{
   value = 0
}
scene_4lvl:add{
  metal_cyl,
  value = 1
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_background = Grid{
   from = {imin-size_pml,jmin-size_pml,kmin-size_pml},
   to = {imax+size_pml,jmax+size_pml,kmax+size_pml}
}
grid_metal_sub = Grid{
   from = {imin-size_pml,jmin-size_pml,1},
   to = {imax+size_pml,jmax+size_pml,kmin-size_pml},
}
--[[
grid_metal_cyl = Grid{
   from = {imin,jmin,-hspacer+2},
   to = {imax,jmax,-cyl_diameter-hspacer-2},
}
--]]
grid_4lvl = Grid{
   from = {-cyl_radius-2,-cyl_length-2,cyl_diameter+hspacer+1},
   to = {cyl_radius+2,cyl_length+2,hspacer-1}
}
grid_xy0 = Grid{
   from = {-cyl_radius-2,0,cyl_diameter+hspacer+1},
   to = {cyl_radius+2,0,hspacer-1}
}
grid_xy1 = Grid{
   from = {-cyl_radius-2,20,cyl_diameter+hspacer+1},
   to = {cyl_radius+2,20,hspacer-1}
}


-- coarse grid for .VTK-preview
grid_prev_background = Grid{
   yee=false,
   from = {imin-size_pml, jmin-size_pml, kmin-size_pml},
   to = {imax+size_pml, jmax+size_pml, kmax+size_pml},
}
grid_prev_metal_sub = Grid{
   yee = false,
   from = {imin-size_pml,jmin-size_pml,1},
   to = {imax+size_pml,jmax+size_pml,kmin-size_pml},
}
grid_prev_4lvl = Grid{
   yee = false,
   from = {-cyl_radius-2,-cyl_length-2,cyl_diameter+hspacer+1},
   to = {cyl_radius+2,cyl_length+2,hspacer-1}
}

cfg:CREATE_GEO{
   "background",
   scene=scene_background,
   grid=grid_background,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_GEO{
   "metal_sub",
   scene=scene_metal_sub,
   grid=grid_metal_sub,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_GEO{
   "metal_cyl",
   scene=scene_metal_cyl,
   grid=grid_metal_cyl,
   method="default",
   comps=3,
   silent=false,
   on=false
}
cfg:CREATE_GEO{
   "4lvl",
   scene=scene_4lvl,
   grid=grid_4lvl,
   method="default",
   comps=3,
   silent=false,
--   on=false
}

cfg:CREATE_GEO{
   "xy0",
   scene=scene_4lvl,
   grid=grid_xy0,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_GEO{
   "xy1",
   scene=scene_4lvl,
   grid=grid_xy1,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}

cfg:CREATE_PREVIEW{
   "background",
   scene=scene_background,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "4lvl",
   scene=scene_4lvl,
   grid=grid_prev_4lvl,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal_sub",
   scene=scene_metal_sub,
   grid=grid_prev_metal_sub,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal_cyl",
   scene=scene_metal_cyl,
   grid=grid_prev_background,
   on=false
}

