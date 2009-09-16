scene_np = Scene{
   value =0.
}
for i,v in ipairs(rnp) do
   sphere1 = Sphere{
      at={inp[i],jnp[i],knp[i]},
      radius=rnp[i]
   }
   scene_cap_inf:add{
      sphere1,
      value = eps_infDL
   }
   scene_np:add{
      sphere1,
      value = 1   -- filling factor for metal!!!!!!!!!!!!
   }
end

table.sort(rnp); table.sort(inp); table.sort(jnp); table.sort(knp)
maxntab = table.maxn(rnp);

grid_np = Grid{
   from={inp[1]-rnp[maxntab]-2,jnp[1]-rnp[maxntab]-2,knp[1]-rnp[maxntab]-2},
   to={inp[maxntab]+rnp[maxntab]+2,jnp[maxntab]+rnp[maxntab]+2,knp[maxntab]+rnp[maxntab]+2}
}
grid_prev_np = Grid{
   yee = false,
   from={inp[1]-rnp[maxntab]-2,jnp[1]-rnp[maxntab]-2,knp[1]-rnp[maxntab]-2},
   to={inp[maxntab]+rnp[maxntab]+2,jnp[maxntab]+rnp[maxntab]+2,knp[maxntab]+rnp[maxntab]+2}
}

cfg:CREATE_GEO{
   "np",
   scene=scene_np,
   grid=grid_np,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "np",
   scene=scene_np,
   grid=grid_prev_np,
   on=true
}
cfg:CREATE_GEO{
   "np_inf",
   scene=scene_cap_inf,
   grid=grid_np,
   method="default",
   comps=3,
   silent=false,
--   on=false
}
cfg:CREATE_PREVIEW{
   "np_inf",
   scene=scene_cap_inf,
   grid=grid_prev_np,
   on=false
}


if (mat == 'gold' or mat == 'silver' or mat == 'carbon') then
   cfg:MAT{                               -- define material with frequency/time-dependent answer
      DRUDE{                              -- define a Drude material
         invlambdapl = conv*real_omegaDL/frequ_factor, -- inverse wavelength of plasma frequency
         gammapl = conv*real_gammaDL/frequ_factor,     -- plasma decay rate
         order = 2
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "np" }
      }
   }
   cfg:MAT{
      LORENTZ{                            -- define a Lorentzian material
         invlambdal = conv*real_omegaL/frequ_factor,   -- inverse plasma wavelength
         gammal = conv*real_gammaL/frequ_factor,       -- resonance width
         deltaepsl = deltaepsl            -- delta epsilon
      },
      REG{                                -- region where material is defined
         LOAD_GEO{ "np" }
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
         LOAD_GEO{ "np" }
      }
   }
end
