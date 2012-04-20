-- basic objects
box_coating = Box{
  from = {imin-1,jmin-1,-hmetal},
  to = {imax+1,jmax+1,-hcoating-hmetal}
}
box_metal1 = Box{
  from = {imin-1,jmin-1,0},
  to = {imax+1,jmax+1,-hmetal}
}
box_spacer = Box{
  from = {imin-1,jmin-1,hspacer},
  to = {imax+1,jmax+1,0}
}
box_metal2 = Box{
  from = {imin-1,jmin-1,hmetal+hspacer},
  to = {imax+1,jmax+1,hspacer}
}
if (hcoating ~= 0 ) then
  box_substrate = Box{
    from = {imin-1,jmin-1,hmetal+hspacer+hcoating},
    to = {imax+1,jmax+1,hmetal+hspacer}
  }
else
  box_substrate = Box{
    from = {imin-1,jmin-1,kmax+size_pml+2},
    to = {imax+1,jmax+1,hmetal+hspacer}
  }
end

-- perforation
box_perf = Box{
  from = {-xperf,-yperf,hmetal+hspacer},
  to = {xperf,yperf,-hcoating-hmetal}
}
if (rounding ~= 0) then
  box_cutout_perf1 = Box{
    from = {xperf-rounding,yperf-rounding,hmetal+hspacer},
    to = {xperf,yperf,-hcoating-hmetal}
  }
  cylinder_perf1 = Cylinder{
    at = {xperf-rounding,yperf-rounding,math.floor(-(hcoating-hspacer)/2+.5)},
    radius = rounding,
    height = 2*hmetal+hcoating+hspacer
  }
  box_cutout_perf2 = Box{
    from = {-xperf+rounding,yperf-rounding,hmetal+hspacer},
    to = {-xperf,yperf,-hcoating-hmetal}
  }
  cylinder_perf2 = Cylinder{
    at = {-xperf+rounding,yperf-rounding,math.floor(-(hcoating-hspacer)/2+.5)},
    radius = rounding,
    height = 2*hmetal+hcoating+hspacer
  }
  box_perf = BinaryAndNot{box_perf,box_cutout_perf1}
  box_perf = BinaryOr{box_perf,cylinder_perf1}
  box_perf = BinaryAndNot{box_perf,box_cutout_perf2}
  box_perf = BinaryOr{box_perf,cylinder_perf2}
  box_perf = BinaryAnd{box_perf,Transform{box_perf,angle=180}}
end
-- complex objects
box_metal1 = BinaryAndNot{box_metal1,box_perf}
box_spacer = BinaryAndNot{box_spacer,box_perf}
box_metal2 = BinaryAndNot{box_metal2,box_perf}
perf = box_perf
if (hcoating ~= 0) then
  perf = BinaryOr{box_perf,box_coating}
end

-- move whole structure so that first interface ist at z=0
box_metal1 = Transform{box_metal1,move={0,0,hshift}}
box_spacer = Transform{box_spacer,move={0,0,hshift}}
box_metal2 = Transform{box_metal2,move={0,0,hshift}}
perf = Transform{perf,move={0,0,hshift}}
box_perf = Transform{box_perf,move={0,0,hshift}}
box_substrate = Transform{box_substrate,move={0,0,hshift}}
if (hcoating ~= 0 ) then
  perf = BinaryOr{perf,box_substrate}
end

spacer = box_spacer
metal2 = box_metal2

scene_background = Scene{
   value = n_bg^2
}
scene_background:add{
   box_metal1,
   value = eps_infDL
}
scene_background:add{
   spacer,
   value = n_spacer^2
}
scene_background:add{
   metal2,
   value = eps_infDL
}
scene_background:add{
   perf,
   value = n_perf^2
}
if hcoating == 0 then
scene_background:add{
   box_substrate,
   value = n_substrate^2
}
end

scene_metal1 = Scene{
   value = 0
}
scene_metal1:add{
   box_metal1,
   value = 1
}
scene_metal2 = Scene{
   value = 0
}
scene_metal2:add{
   metal2,
   value = 1
}
scene_4lvl = Scene{
   value = 0
}
if active == 1 then
   scene_4lvl:add{
      perf,
      value = 1
   }
   scene_4lvl:add{
      spacer,
      value = 1
   }
elseif active == 2 then
   scene_4lvl:add{
      perf,
      value = 1
   }
elseif active == 3 then
   scene_4lvl:add{
      spacer,
      value = 1
   }
else
   print("ERROR: no active material")
end

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_background = Grid{
   from = {imin,jmin,kmax+size_pml},
   to = {imax,jmax,-2}
}
grid_metal1 = Grid{
   from = {imin,jmin,hcoating-2},
   to = {imax,jmax,hcoating+hmetal+2},
}
grid_metal2 = Grid{
   from = {imin,jmin,hcoating+hmetal+hspacer-2},
   to = {imax,jmax,hcoating+hmetal+hshift2+2},
}
grid_4lvl = Grid{
   from = {imin,jmin,-2},
   to = {imax,jmax,2*hcoating+hmetal+hshift2+2}
}
-- coarse grid for .VTK-preview
grid_prev_background = Grid{
   yee=false,
   from = {imin, jmin, kmin-size_pml},
   to = {imax, jmax, kmax+size_pml},
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
   "metal1",
   scene=scene_metal1,
   grid=grid_metal1,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_GEO{
   "metal2",
   scene=scene_metal2,
   grid=grid_metal2,
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
   "4lvl",
   scene=scene_4lvl,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal1",
   scene=scene_metal1,
   grid=grid_prev_background,
--   on=false
}
cfg:CREATE_PREVIEW{
   "metal2",
   scene=scene_metal2,
   grid=grid_prev_background,
--   on=false
}
