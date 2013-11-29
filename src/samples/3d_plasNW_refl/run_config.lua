--cfg = CONFIG{scenes=false}
cfg = CONFIG{scenes=true}
cyl_rad = select(1,...)
refl = (select(2,...) or 0)
mode = (select(3,...) or 0)
neff_guess = select(4,...)
inj_pos = select(5,...)
end_pos = 30

dofile("scale.lua")

imin = -hdist_tfsf_i
imax = hdist_tfsf_i
jmin = -hdist_tfsf_jn
jmax = hdist_tfsf_jp
kmin = -hdist_tfsf_kn
kmax = hdist_tfsf_kp

dist_end = 10
if (refl=='1') then
  inj_pos = inj_pos - 30
--  inj_pos = -4--kmax-8-dist_end
  box_xy0 = { imin, imax, 1, jmin, jmax, 1, end_pos+30, end_pos+30, 1 }
  box_xy1 = { imin, imax, 1, jmin, jmax, 1, end_pos+60, end_pos+60, 1 }
  box_xy2 = { imin, imax, 1, jmin, jmax, 1, 0, 0, 1 }
  box_xy3 = { imin, imax, 1, jmin, jmax, 1, math.floor(kmin/2), math.floor(kmin/2), 1 }
  box_xy4 = { imin, imax, 1, jmin, jmax, 1, kmin+2, kmin+2, 1 }
else
  inj_pos = kmin+4+dist_end
  box_xy0 = { imin, imax, 1, jmin, jmax, 1, inj_pos-4, inj_pos-4, 1 }
  box_xy1 = { imin, imax, 1, jmin, jmax, 1, inj_pos+4, inj_pos+4, 1 }
  box_xy2 = { imin, imax, 1, jmin, jmax, 1, 0, 0, 1 }
  box_xy3 = { imin, imax, 1, jmin, jmax, 1, math.floor(kmax/2), math.floor(kmax/2), 1 }
  box_xy4 = { imin, imax, 1, jmin, jmax, 1, kmax-2, kmax-2, 1 }
end
print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 3,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = {0,ncycles},                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin-size_pml, imax+size_pml },   -- range of computational window in x direction
   jrange = { jmin-size_pml, jmax+size_pml },   -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

}

cfg:CHECKPOINT{

   load = false,
   save = false,
   detail = 3

}

--- CREATE SCENE
dofile("run_geo.lua")

--- fire up matlab mode calculator!
nrefr = (neff_guess or n_refr)
if (mode=='1') then
  print("> luacfg mode.lua "..tostring(pinv_wl).." "..tostring(nrefr).." 1")
  silent=1
  --os.execute("~/fdtd3d-par/trunk/src/build-gcc/luacfg mode.lua "..tostring(pinv_wl).." "..tostring(nrefr).." 1")
  dofile("mode.lua")
end
dofile("neff.lua")

--- FDTD Definition
cfg:FDTD{

   EPSILON{
      REG{
         BOX( box_sub ),
         LOAD_GEO{ "background" },
      },
      on = true
   },

   OUT{
      file = { "VTK", "epsxz" },
      type = { "Eps", "N" },
      time = { 0, 0, 1 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      },
      on = false
   },
   OUT{  
      file = { "GPL", "en_sum" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            { imin, imax, 1, 0, cyl_diameter+hspacer, 1, jmin, jmax, 1 }
         } 
      }
   },
   OUT{
      file = { "GPL", "en_xy0" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy0" },
      }
   },
   OUT{
      file = { "GPL", "en_xy1" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy1" },
      }
   },
   OUT{
      file = { "GPL", "en_xy2" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy2" },
      }
   },
   OUT{
      file = { "GPL", "en_xy3" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy3" },
      }
   },
   OUT{
      file = { "GPL", "en_xy4" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy4" },
      }
   },
   OUT{
      file = { "GPL", "Sz0" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy0
         }
      }
   },
   OUT{
      file = { "GPL", "Sz1" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy1
         }
      }
   },
   OUT{
      file = { "GPL", "Sz2" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy2
         }
      }
   },
   OUT{
      file = { "GPL", "Sz3" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy3
         }
      }
   },
   OUT{
      file = { "GPL", "Sz4" },
      type = { "Sz", "S", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy4
         }
      }
   },

   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { ncycles, ncycles, 1 },
      REG{ 
         BOX{
            { imin, imax, 1, jmin, jmax, 1, kmin, kmax, 1}
         }
      },    
      on = false
   },
   OUT{
      file = { "VTK", "E_xy_inj" },
      type = { "E", "N" },
      time = { 1550, 1650, 1 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, 0, 1 }
         }
      },
      on = injection_vtk_on
   },
   OUT{
      file = { "VTK", "H_xy_inj" },
      type = { "H", "N" },
      time = { 1550, 1650, 1 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, 0, 1 }
         }
      },
      on = injection_vtk_on
   },
   OUT{
      file = { "VTK", "E_xy0" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            box_xy0
         }
      },
      on = false--field_vtk_on
   },
   OUT{
      file = { "VTK", "E_xy1" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            box_xy1
         }
      },
      on = false--field_vtk_on
   },
   OUT{
      file = { "VTK", "E_xy2" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            box_xy2
         }
      },
      on = false--field_vtk_on
   },
   OUT{
      file = { "VTK", "E_xy3" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            box_xy3
         }
      },
      on = false--field_vtk_on
   },
   OUT{
      file = { "VTK", "E_xy4" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            box_xy4
         }
      },
      on = false--field_vtk_on
   },
   OUT{
      file = { "VTK", "E_yz" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      },
      on = field_vtk_on
   },
   OUT{
      file = { "VTK", "E_xz" },
      type = { "E", "N" },
      time = { 500, ncycles, 500 },
      REG{
         BOX{
         { imin, imax, step_dft, 0, 0, 1, kmin, kmax, step_dft }
         }
      },
      on = field_vtk_on
   },

   OUT{
      file = { "VTK", "E_xz" },
      type = { "E", "N" },
      time = { ncyc_probe_start-3000, ncyc_probe_start+3000, 100 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "E_yz" },
      type = { "E", "N" },
      time = { ncyc_probe_start-3000,ncyc_probe_start+3000, 100 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      },
      on = false
   },

 
}

--- BOUND Definition

cfg:BOUND{

   config = { 4, 4, 4, 4, 4, 4 },

   CPML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
      alpha = 0.05,
      alphapot = 3.2
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFBOX{
      invlambda = inv_wl,
      amplitude = ampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(widthl/inv_wl/dt+.5),
         offset=math.floor(offsetl/inv_wl/dt+.5),
         attack=math.floor(attackl/inv_wl/dt+.5),
         sustain=math.floor((sustainl+0*10000)/inv_wl/dt+.5),
         decay=math.floor(decayl/inv_wl/dt+.5)
      },
      planewave = { phi=-90, theta=90.0, psi=psi, nrefr=nrefr },
      config = {0,0,0,1,0,0}
   },
   REG{
      BOX{
         { imin-size_pml-1, imax+size_pml+1, 1, jmin, jmax, 1, kmin-size_pml-1, kmax+size_pml+1, 1 }
      }
   },
   on = pump_on
}

cfg:SRC{
   TFSFBOX{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      planewave = { phi=-90, theta=90, psi=ppsi, nrefr=nrefr },
      config = {0,0,0,1,0,0}
   },
   REG{
      BOX{
         { imin-size_pml-1, imax+size_pml+1, 1, jmin, jmax, 1, kmin-size_pml-1, kmax+size_pml+1, 1 }
      }
   },
   on = probe_on
}
cfg:SRC{
   TFSFINJ{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      planewave = { phi=0, theta=0, psi=ppsi-180, nrefr=nrefr },
--      config = {1,1,1,1,1,0}
   },
   REG{
      LOAD{ "tfsfex_re" }
   },
   on = probe_on_injx,
}
   cfg:SRC{
   TFSFINJ{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      alpha = -90,
      planewave = { phi=0, theta=0, psi=ppsi-180, nrefr=nrefr },
--      config = {1,1,1,1,1,0}
   },
   REG{
      LOAD{ "tfsfex_im" }
   },
   on = probe_on_injx,
}
cfg:SRC{
   TFSFINJ{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      planewave = { phi=0, theta=0, psi=ppsi-270, nrefr=nrefr },
--      config = {1,1,1,1,1,0}
   },
   REG{
      LOAD{ "tfsfey_re" }
   },
   on = probe_on_injy
}
cfg:SRC{
   TFSFINJ{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      alpha = -90,
      planewave = { phi=0, theta=0, psi=ppsi-270, nrefr=nrefr },
--      config = {1,1,1,1,1,0}
   },
   REG{
      LOAD{ "tfsfey_im" }
   },
   on = probe_on_injy
}

--- MAT Definition(s)
if (mat == 'gold' or mat == 'silver' or mat == 'silver2' ) then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_sub" }
      }
   }
   if ( mat ~= 'silver2') then
      cfg:MAT{
         LORENTZ{                            -- define a Lorentzian material
            invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
            gammal = conv*real_gammaL/frequ_factor,       -- resonance width
            deltaepsl = deltaepsl            -- delta epsilon
         },
         REG{                                -- region where material is defined
            LOAD_GEO{ "metal_sub" }
         }
      }
   end
end
if (mat == 'silver') then
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL2/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL2/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl2            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_sub" }
      }
   }
end

cfg:MAT{
   LVFOURLVLEP{
      invlambdal = {invlambdala,invlambdalb},
      gammal = {gammala,gammalb},
      dipole12 = dipolea,
      dipole03 = dipoleb,
      rp = pump_rate,
      dens = density,
      start = {0,init_inv,0},
      gamma = {gamma30,gamma32,gamma21,gamma10},
      LFE = localfield,
      volfac = volfactor,
      linefac = linefactor,
   },
   REG{
      LOAD_GEO{ "4lvl" },
   },
   OUT{
      file = { "GPL", "dens_sum" },
      type = { "N", "S", ".F." },
      time = { 0, ncycles, 100 },
      REG{
         BOX{
            { imin, imax, 1, -1, cyl_diameter+hspacer+1, 1, kmin, kmax, 1 }
         }
      },
   },
   OUT{
      file = { "VTK", "dens_xz" },
      type = { "N", "N" },
      time = { 0, ncycles, 1000 },
      REG{
         BOX{
            { imin, imax, 1, -1, cyl_diameter+hspacer+1, 1, 0, 0, 1 }
         }
      },
--      on = false
   },

OUT{
      file = { "VTK", "dens_yz" },
      type = { "N", "N" },
      time = { 0, ncycles, 1000 },
      REG{
         BOX{
            { 0,  0, 1, -1, cyl_diameter+hspacer+1, 1, kmin, kmax, 1 }
         }
      },
--      on = false
   },

   on = on_4lvl
}

-- Energy balance module
if (diagebal) then
  dofile("diagebal.lua")
end

-- Spectrum diagnosis module
if (diagpspec) then
  dofile("diagpspec.lua")
end

-- diagnostics DFT/FFT
if (diagmode) then
  dofile("diagmode.lua")
end

--- CREATE: config.<part>.in
cfg:CREATE()

