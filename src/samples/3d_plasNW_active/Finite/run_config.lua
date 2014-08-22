--cfg = CONFIG{scenes=false}
cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = -hdist_ntff_i
imax = hdist_ntff_i-1
jmin = -hdist_ntff_j
jmax = hdist_ntff_j-1
kmin = -hdist_ntff_kn
kmax = hdist_ntff_kp

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
   dx = { conv*1e-9/1e40, 1, 1, 1 }

}

cfg:CHECKPOINT{

   load = false,
   save = true,
   detail = 3

}

--- CREATE SCENE
dofile("run_geo.lua")

--- FDTD Definition
cfg:FDTD{

   EPSILON{
      REG{
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
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, cyl_diameter+hspacer, 1 }
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
      file = { "VTK", "E_xz_late" },
      type = { "E", "N" },
      time = { 50000-400, 50000, 10 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      },
--      on = false
   },
   OUT{
      file = { "VTK", "E_xz2_late" },
      type = { "E", "N" },
      time = { 50000-400, 50000, 10 },
      REG{
         BOX{
            { imin, imax, 1, 6, 6, 1, kmin, kmax, 1 }
         }
      },
--      on = false
   },
   OUT{
      file = { "VTK", "E_xz3_late" },
      type = { "E", "N" },
      time = { 50000-400, 50000, 10 },
      REG{
         BOX{
            { imin, imax, 1, 77, 77, 1, kmin, kmax, 1 }
         }
      },
--      on = false
   },
   OUT{
      file = { "VTK", "E_yz_late" },
      type = { "E", "N" },
      time = { 50000-400, 50000, 10 },
      REG{
         BOX{
            { 0, 0, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      },
--      on = false
   },
   OUT{
      file = { "VTK", "E_xy_late" },
      type = { "E", "N" },
      time = { 50000-400, 50000, 10 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 1, 1, 1 }
         }
      },
--      on = false
   },

   OUT{
      file = { "VTK", "E_xz" },
      type = { "E", "N" },
      time = { 20, 2000, 20 },
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
      time = { 20, 2000, 20 },
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
      amplitude = ampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)*1e20),
      pulse = {
         shape="Sech",
         width=math.floor(widthl*resolution/dt+.5),
         offset=math.floor(offsetl*resolution/dt+.5),
         attack=math.floor(attackl*resolution/dt+.5),
         sustain=math.floor((sustainl+0*10000)*resolution/dt+.5),
         decay=math.floor(decayl*resolution/dt+.5)
      },
      planewave = { phi=0, theta=0.0, psi=psi, nrefr=nrefr },
      config = {1,1,1,1,1,0}
   },
   REG{
      BOX{
         {-hdist_tfsf_i-1,hdist_tfsf_i+1,1,-hdist_tfsf_j-1,hdist_tfsf_j+1,1,-hdist_tfsf_kn,10,1}
      }
   },
   on = pump_on
}

cfg:SRC{
   TFSFBOX{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)*1e20),
      pulse = {
         shape="Sech",
         width=math.floor(pwidthl*resolution/dt+.5),
         offset=math.floor(poffsetl*resolution/dt+.5),
         attack=math.floor(pattackl*resolution/dt+.5),
         sustain=math.floor(psustainl*resolution/dt+.5),
         decay=math.floor(pdecayl*resolution/dt+.5)
      },
      planewave = { phi=0, theta=0, psi=ppsi, nrefr=nrefr },
      config = {1,1,1,1,1,0}
   },
   REG{
      BOX{
         {-hdist_tfsf_i-1,hdist_tfsf_i+1,1,-hdist_tfsf_j-1,hdist_tfsf_j+1,1,-hdist_tfsf_kn,10,1}
      }
   },
   on = probe_on
}


--- MAT Definition(s)
--[[
if (mat == 'gold' or mat == 'silver' ) then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_cyl" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_cyl" }
      }
   }
end
if (mat == 'silver') then
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL2/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL2/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl2            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_cyl" }
      }
   }
end
--]]
if (mat == 'gold' or mat == 'silver') then
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
            { imin, imax, 1, jmin, jmax, 1, -1, cyl_diameter+hspacer+1, 1 }
         }
      },
   },
   OUT{
      file = { "VTK", "dens" },
      type = { "dens", "N" },
      time = { ncycles, ncycles, 1 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, -1, cyl_diameter+hspacer+1, 1 }
         }
      },
      on = false
   },
--   on = false
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

