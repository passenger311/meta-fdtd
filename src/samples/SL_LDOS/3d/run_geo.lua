-- basic objects
box_metal_sub = Box{
  from = {imin-1,jmin-1,0},
  to = {imax+1,jmax+1,hmetal_sub+2}
}
spacer = Box{
  from = {imin-1,jmin-1,-hspacer},
  to = {imax+1,jmax+1,0}
}
box_spacer = Box{
  from = {-hdist_tfsf_i,-hdist_tfsf_j,-hspacer},
  to = {hdist_tfsf_i,hdist_tfsf_j,0}
}
box_metal_bar = Box{
  from = {imin-1,jmin-1,-bar_height-hspacer},
  to = {imax+1,jmax+1,-hspacer}
}
box_top = Box{
  from = {imin-1,jmin-1,-bar_height-hspacer},
  to = {imax+1,jmax+1,-bar_height-hspacer-top}
}

-- complex objects

metal_bar = box_metal_bar

scene_background = Scene{
   value = n_bg^2
}
scene_background:add{
   box_metal_sub,
   value = eps_infDL
}
scene_background:add{
   spacer,
   value = n_spacer^2
}
scene_background:add{
   metal_bar,
   value = eps_infDL
}
scene_background:add{
   box_top,
   value = n_spacer^2
}
scene_metal_sub = Scene{
   value = 0
}
scene_metal_sub:add{
   box_metal_sub,
   value = 1
}
scene_metal_bar = Scene{
   value = 0
}
scene_metal_bar:add{
   metal_bar,
   value = 1
}
scene_4lvl = Scene{
   value = 0
}
scene_4lvl:add{
  spacer,
  value = 1
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_background = Grid{
   from = {imin,jmin,-bar_height-hspacer-top-2},
   to = {imax,jmax,hmetal_sub}
}
grid_metal_sub = Grid{
   from = {imin,jmin,-2},
   to = {imax,jmax,hmetal_sub+2},
}
grid_metal_bar = Grid{
   from = {imin,jmin,-hspacer+2},
   to = {imax,jmax,-bar_height-hspacer-2},
}
grid_4lvl = Grid{
   from = {imin,jmin,2},
   to = {imax,jmax,-hspacer-2}
}
-- coarse grid for .VTK-preview
grid_prev_background = Grid{
   yee=false,
   from = {0, jmin, kmin-size_pml},
   to = {0, jmax, kmax},
}
grid_prev_background2 = Grid{
   yee=false,
   from = {imin, 0, kmin-size_pml},
   to = {imax, 0, kmax},
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
   "metal_bar",
   scene=scene_metal_bar,
   grid=grid_metal_bar,
   method="default",
   comps=3,
   silent=false,
--   on=false
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
cfg:CREATE_PREVIEW{
   "background",
   scene=scene_background,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "background2",
   scene=scene_background,
   grid=grid_prev_background2,
--   on=false
}
cfg:CREATE_PREVIEW{
   "4lvl",
   scene=scene_4lvl,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal_sub",
   scene=scene_metal_sub,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal_bar",
   scene=scene_metal_bar,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal_bar2",
   scene=scene_metal_bar,
   grid=grid_prev_background2,
--   on=false
}

print(kmax,kmin-size_pml)
print(hmetal_sub+2,0,-hspacer,-bar_height-hspacer,-bar_height-hspacer-top)

