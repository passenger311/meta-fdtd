-- basic objects
box_metal_sub = Box{
  from = {imin-size_pml-1,0,kmin-size_pml-1},
  to = {imax+size_pml+size_pml+1,-hmetal_sub-size_pml-1,kmax+size_pml+1}
}
spacer = Box{
  from = {imin-size_pml-1,0,kmin-size_pml-1},
  to = {imax+size_pml+1,hspacer,kmax+size_pml+1}
}
box_sub = {}
if (mat~='dielectric') then
  table.insert( box_sub, {imin-size_pml, imax+size_pml, 1, jmin-size_pml, 0, 1, kmin-size_pml, kmax+size_pml, 1,":",eps_infDL,eps_infDL,eps_infDL })
end
-- complex objects

ZnO_cyl = Cylinder{
  radius = cyl_radius1,
  height = mult_factor*2*cyl_length,
  at = {0, 0, 0}
}
--ZnO_cyl = Transform{ZnO_cyl,angle=90,axis={1, 0, 0}}
ZnO_cyl = Transform{ZnO_cyl,move={0,cyl_radius1+hspacer,-10*cyl_length+inj_pos+2}} --delete last bit for active

scene_bg = Scene{
   value = n_bg^2
}
scene_bg:add{
   box_metal_sub,
   value = eps_infDL
}
scene_bg:add{
   spacer,
   value = n_spacer^2
}
scene_bg:add{
   ZnO_cyl,
   value = n_nanowire^2
}

scene_metal_sub = Scene{
   value = 0
}
scene_metal_sub:add{
   box_metal_sub,
   value = 1
}
scene_4lvl = Scene{
   value = 0
}
scene_4lvl:add{
  ZnO_cyl,
  value = 1
}

-- specify a grid for the scene (only objects that are inside the grid will be part of the geometry)
grid_background = Grid{
   from = {imin-size_pml,-2,kmin-size_pml},
   to = {imax+size_pml,jmax+size_pml,kmax+size_pml}
}
grid_metal_sub = Grid{
   from = {imin-size_pml,1,kmin-size_pml},
   to = {imax+size_pml,jmin-size_pml,kmax+size_pml},
}
grid_4lvl = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,-cyl_length-2},
   to = {cyl_radius+2,hspacer-1,cyl_length+2}
}
grid_xy0 = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,inj_pos-2},
   to = {cyl_radius+2,hspacer-1,inj_pos-2}
}
grid_xy1 = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,math.floor(inj_pos/2)},
   to = {cyl_radius+2,hspacer-1,math.floor(inj_pos/2)}
}
grid_xy2 = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,-2},
   to = {cyl_radius+2,hspacer-1,-2}
}
grid_xy3 = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,-math.floor(inj_pos/2)},
   to = {cyl_radius+2,hspacer-1,-math.floor(inj_pos/2)}
}
grid_xy4 = Grid{
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,-inj_pos+2},
   to = {cyl_radius+2,hspacer-1,-inj_pos+2}
}

-- coarse grid for .VTK-preview
grid_prev_background = Grid{
   yee=false,
   from = {imin-size_pml, -2, kmin-size_pml},
   to = {imax+size_pml, jmax+size_pml, kmax+size_pml},
}
grid_prev_metal_sub = Grid{
   yee = false,
   from = {imin-size_pml,1,kmin-size_pml},
   to = {imax+size_pml,jmin-size_pml,kmax+size_pml},
}
grid_prev_4lvl = Grid{
   yee = false,
   from = {-cyl_radius-2,cyl_diameter+hspacer+1,-cyl_length-2},
   to = {cyl_radius+2,hspacer-1,cyl_length+2}
}

cfg:CREATE_GEO{
   "background",
   scene=scene_bg,
   grid=grid_background,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_GEO{
   "metal_sub",
   scene=scene_metal_sub,
   grid=grid_metal_sub,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_GEO{
   "ZnO_cyl",
   scene=scene_ZnO_cyl,
   grid=grid_ZnO_cyl,
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
   on=geo_on
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
cfg:CREATE_GEO{
   "xy2",
   scene=scene_4lvl,
   grid=grid_xy2,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_GEO{
   "xy3",
   scene=scene_4lvl,
   grid=grid_xy3,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_GEO{
   "xy4",
   scene=scene_4lvl,
   grid=grid_xy4,
   method="default",
   comps=3,
   silent=false,
   on=geo_on
}
cfg:CREATE_PREVIEW{
   "background",
   scene=scene_bg,
   grid=grid_prev_background,
   on=prev_on
}
cfg:CREATE_PREVIEW{
   "4lvl",
   scene=scene_4lvl,
   grid=grid_prev_4lvl,
   on=prev_on
}
cfg:CREATE_PREVIEW{
   "metal_sub",
   scene=scene_metal_sub,
   grid=grid_prev_metal_sub,
   on=prev_on
}
cfg:CREATE_PREVIEW{
   "ZnO_cyl",
   scene=scene_ZnO_cyl,
   grid=grid_prev_background,
   on=false
}


--Injection plane calculation
d_om=0.001
if ( mat == 'silver') then
  om = 2.9979e5/real_wavelength
  e_inf = eps_infDL
  oDL = real_omegaDL
  gDL = real_gammaDL/2/pi
  oL1 = real_omegaL
  gL1 = real_gammaL/2/pi
  dL1 = deltaepsl
  oL2 = real_omegaL2
  gL2 = real_gammaL2/2/pi
  dL2 = deltaepsl2
  eps_metal = e_inf - oDL^2/(om^2+gDL^2) - dL1*oL1^2*(om^2-oL1^2)/((om^2-oL1^2)^2+4*gL1^2*om^2) - dL2*oL2^2*(om^2-oL2^2)/((om^2-oL2^2)^2+4*gL2^2*om^2)
  eps_metal_im = oDL^2/(om^2+gDL^2)*gDL/om + dL1*oL1^2*(2*gL1*om)/((om^2-oL1^2)^2+4*gL1^2*om^2) + dL2*oL2^2*(2*gL2*om)/((om^2-oL2^2)^2+4*gL2^2*om^2)
  om = om+d_om
  eps_metal2_re = e_inf - oDL^2/(om^2+gDL^2) - dL1*oL1^2*(om^2-oL1^2)/((om^2-oL1^2)^2+4*gL1^2*om^2) - dL2*oL2^2*(om^2-oL2^2)/((om^2-oL2^2)^2+4*gL2^2*om^2)
  deps_metal = eps_metal2_re - eps_metal
  om = om-d_om
  d_om_eps = eps_metal+om*deps_metal/d_om
  print("metal permittivity", eps_metal,eps_metal_im)
  print("d(w eps)/dw", d_om_eps)

end
if ( mat == 'gold') then
  om = 2.9979e5/real_wavelength
  e_inf = eps_infDL
  oDL = real_omegaDL
  gDL = real_gammaDL/2/pi
  oL1 = real_omegaL
  gL1 = real_gammaL/2/pi
  dL1 = deltaepsl
  eps_metal = e_inf - oDL^2/(om^2+gDL^2) - dL1*oL1^2*(om^2-oL1^2)/((om^2-oL1^2)^2+4*gL1^2*om^2)
  eps_metal_im = oDL^2/(om^2+gDL^2)*gDL/om + dL1*oL1^2*(2*gL1*om)/((om^2-oL1^2)^2+4*gL1^2*om^2)
  om = om+d_om
  eps_metal2_re = e_inf - oDL^2/(om^2+gDL^2) - dL1*oL1^2*(om^2-oL1^2)/((om^2-oL1^2)^2+4*gL1^2*om^2)
  deps_metal = eps_metal2_re - eps_metal
  om = om-d_om
  d_om_eps = eps_metal+om*deps_metal/d_om
  print("metal permittivity", eps_metal,eps_metal_im)
  print("d(w eps)/dw", d_om_eps)
end
if ( mat == 'silver2') then
  om = 2.9979e5/real_wavelength
  e_inf = eps_infDL
  oDL = real_omegaDL
  gDL = real_gammaDL/2/pi
  eps_metal = e_inf - oDL^2/(om^2+gDL^2)
  eps_metal_im = oDL^2/(om^2+gDL^2)*gDL/om
  om = om+d_om
  eps_metal2_re = e_inf - oDL^2/(om^2+gDL^2) 
  deps_metal = eps_metal2_re - eps_metal
  om = om-d_om
  d_om_eps = eps_metal+om*deps_metal/d_om
  print("metal permittivity", eps_metal,eps_metal_im)
  print("d(w eps)/dw", d_om_eps)
end
if ( mat == 'dielectric') then
  om = 2.9979e5/real_wavelength
  e_inf = eps_infDL
  eps_metal = e_inf
  eps_metal_im = 0

end
scene_injection_re = Scene{
   value = n_bg^2
}
scene_injection_re:add{
   box_metal_sub,
   value = eps_metal
}
scene_injection_re:add{
   spacer,
   value = n_spacer^2
}
scene_injection_re:add{
   ZnO_cyl,
   value = n_nanowire^2
}
scene_injection_im = Scene{
   value = 0
}scene_injection_im_gain = Scene{
   value = 0
}
if (mat ~= 'dielelectric') then
  scene_injection_im:add{
     box_metal_sub,
     value = eps_metal_im
  }
  scene_injection_im_gain:add{
     ZnO_cyl,
     value = -0.00
  }
end
scene_nanowire = Scene{
   value = 0
}
scene_nanowire:add{
   ZnO_cyl,
   value = n_nanowire^2
}
scene_energy = Scene{
   value = n_bg^2
}
scene_energy:add{
   box_metal_sub,
   value = d_om_eps
}
scene_energy:add{
   spacer,
   value = n_spacer^2
}
scene_energy:add{
   ZnO_cyl,
   value = n_nanowire^2
}


grid_injection = Grid{
   from = {imin-size_pml,jmin-size_pml,inj_pos},
   to = {imax+size_pml,jmax+size_pml,inj_pos},
}
grid_prev_injection = Grid{
   from = {imin-size_pml,jmin-size_pml,inj_pos},
   to = {imax+size_pml,jmax+size_pml,inj_pos},
   yee=false
}
cfg:CREATE_GEO{
   "injection_re",
   scene=scene_injection_re,
   grid=grid_injection,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "injection_re",
   scene=scene_injection_re,
   grid=grid_prev_injection,
   on=false
}
cfg:CREATE_GEO{
   "injection_im_metal",
   scene=scene_injection_im,
   grid=grid_injection,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_GEO{
   "injection_im_gain",
   scene=scene_injection_im_gain,
   grid=grid_injection,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "injection_im",
   scene=scene_injection_im,
   grid=grid_prev_injection,
   on=false
}
cfg:CREATE_GEO{
   "nanowire",
   scene=scene_nanowire,
   grid=grid_injection,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "nanowire",
   scene=scene_nanowire,
   grid=grid_prev_injection,
   on=false
}
cfg:CREATE_GEO{
   "energy",
   scene=scene_energy,
   grid=grid_injection,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "energy",
   scene=scene_energy,
   grid=grid_prev_injection,
   on=false
}
