--cfg = CONFIG{scenes=false}
cfg = CONFIG{scenes=true}
taskid = select(1,...)
dofile("scale.lua")

imin = -hdist_ntff_i
imax = hdist_ntff_i-1
jmin = -hdist_ntff_j
jmax = hdist_ntff_j-1
kmin = -hdist_ntff_kn-size_pad
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
   irange = { imin, imax },   -- range of computational window in x direction
   jrange = { jmin, jmax },   -- -"- in y direction
   krange = { kmin-size_pml, kmax },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

}

cfg:CHECKPOINT{

   load = false,
   save = false,
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
      file = { "GPL", "E" },
      type = { "E", "S", ".F." },
      time = { 0, ncycles, 1 },
      REG{
         POINT{
            { 0, 0, -hspacer+1+taskid }
         }
      },
      on = true
   },
   OUT{  
      file = { "GPL", "en_sum" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, bar_height+hspacer, 1 }
         } 
      },
      on = false
   },

   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { ncycles, ncycles, 1 },
      REG{ 
         BOX{
            { imin, imax, 1, jmin, jmax, 1, -hspacer, bar_height+2*hspacer, 1 }
         }
      },    
      on = false
   },

   OUT{
      file = { "VTK", "E_xz_late" },
      type = { "E", "N" },
      time = { ncycles-400, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 1, 1, 1, 1, kmin-size_pml, kmax+size_pml, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "E_yz_late" },
      type = { "E", "N" },
      time = { ncycles-400, ncycles, 10 },
      REG{
         BOX{
            { 1, 1, 1, jmin, jmax, 1, kmin-size_pml, kmax+size_pml, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "E_xz" },
      type = { "E", "N" },
      time = { math.floor(pattackl*resolution/dt+.5)-4*math.floor(pwidthl*resolution/dt+.5), math.floor(pattackl*resolution/dt+.5)+4*math.floor(pwidthl*resolution/dt+.5), 20 },
      REG{
         BOX{
            { imin, imax, 1, 1, 1, 1, kmin-size_pml, kmax, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "E_yz" },
      type = { "E", "N" },
      time = { math.floor(pattackl*resolution/dt+.5)-4*math.floor(pwidthl*resolution/dt+.5), math.floor(pattackl*resolution/dt+.5)+4*math.floor(pwidthl*resolution/dt+.5), 20 },
      REG{
         BOX{
            { 1, 1, 1, jmin, jmax, 1, kmin-size_pml, kmax, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "E_xy" },
      type = { "E", "N" },
      time = { math.floor(pattackl*resolution/dt+.5)-4*math.floor(pwidthl*resolution/dt+.5), math.floor(pattackl*resolution/dt+.5)+4*math.floor(pwidthl*resolution/dt+.5), 20 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 1, 1, 1 }
         }
      },
      on = false
   },
 
}

--- BOUND Definition

cfg:BOUND{

   config = { 4, 4, 4, 4, 4, 0 },

   CPML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
      alpha = 0.05,
      alphapot = 3
   }
}
--[[
cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 0 },

   PML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
--      alpha = 0.05,
--      alphapot = 1
   }
}
--]]
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
   HARDJ{
      invlambda = pinv_wl,
--      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      amplitude = 1e6,
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl*resolution/dt+.5),
         offset=math.floor(poffsetl*resolution/dt+.5),
         attack=math.floor(pattackl*resolution/dt+.5),
         sustain=math.floor(psustainl*resolution/dt+.5),
         decay=math.floor(pdecayl*resolution/dt+.5)
      },
--      planewave = { on=false, phi=0, theta=0, psi=ppsi, nrefr=nrefr },
--      planewave.on={false},
--      config = {0,0,0,0,1,1}
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -hspacer+1+taskid, -hspacer+1+taskid, 1, ":", 1, 0, 0}
      },
--      BOX{
--         {-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf,hdist_tfsf,1}
--      }
   },
   on = probe_on
}


--- MAT Definition(s)

if (mat == 'gold' or mat == 'silver' or mat =='ITO' ) then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_bar" }
      }
   }
end
if (mat == 'gold' or mat == 'silver' ) then
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "metal_bar" }
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
         LOAD_GEO{ "metal_bar" }
      }
   }
end
if (mat == 'gold' or mat == 'silver' or mat == 'ITO') then
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
end
if (mat == 'gold' or mat == 'silver' ) then
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
   LVFOURLVL{
      invlambdal = {invlambdala,invlambdalb},
      gammal = {gammala,gammalb},
      dipole12 = dipolea,
      dipole03 = dipoleb,
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
            { imin, imax, 1, jmin, jmax, 1, -1, bar_height+hspacer+1, 1 }
         }
      },
   },
   OUT{
      file = { "VTK", "dens" },
      type = { "dens", "N" },
      time = { ncycles, ncycles, 1 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, -1, bar_height+hspacer+1, 1 }
         }
      },
      on = false
   },
   on = false
}


--- CREATE: config.<part>.in
cfg:CREATE()

