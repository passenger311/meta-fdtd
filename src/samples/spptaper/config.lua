cfg = CONFIG{scenes=true}    -- scences takes either true or false 

TASKID = select(1,...)

dofile("scale.lua")          -- geometrical and material parameters are defined in scale file

-- grid dimensions: i,j,k

imax = round( hwidth )   -- i -> width
jmax = round( height )   -- j -> height
kmax = round( length )   -- k -> direction of propagation

--- inner computational domain

imin = 0
jmin = 0
kmin = 0

print("irange (grid) = ", imin, imax)
print("jrange (grid) = ", jmin, jmax)
print("krange (grid) = ", kmin, kmax)
print("ncyc          = ", ncyc)

print("irange (real) = ", imin*real_dx, imax*real_dx)
print("jrange (real) = ", jmin*real_dx, jmax*real_dx)
print("krange (real) = ", kmin*real_dx, kmax*real_dx)

--- fix some parameters

pad1=0

imin0 = imin
imax0 = imax +pad1 + cpml
jmin0 = jmin -pad1 - cpml
jmax0 = jmax +pad1 + cpml
kmin0 = kmin -pad1 - cpml
kmax0 = kmax +pad1 + cpml


--- dielectric structure!

sc_eps = Scene{ value=eps_bg }  

-- ... sio2 base

box_bsio2 = Box{
   from={-imax0-100,jmin0-100,kmin0-100},
   to={imax0+100,height_bsio2,kmax0+100}
}

sc_eps:add{ box_bsio2, depth=1, value=eps_sio2 }

-- ... si base

box_bsi = Box{ 
   from={-imax0-100,height_bsio2,kmin0-100},
   to={imax0+100,height_bsio2+height_bsi,kmax0+100}
}

sc_eps:add{ box_bsi, depth=1, value=eps_si }     

-- .. wg structure

height1 = height_bsio2+height_bsi
height2 = height_bsio2+height_bsi+height_wg

l = length_it / math.cos(alpha)
w = l * math.tan(alpha)


box_is1 = Transform{ 
   Box{
      from = { 0  , height1,  0 },
      to = { w, height2, l }
   },
   axis = { 0, 1, 0},
   angle = - alpha  * 180 / math.pi,
   move = { hwidth_i, 0, length_i  }
}

box_is2 = Transform{ 
   Box{
      from = { -w  , height1,  0 },
      to = { 0 , height2, l }
   },
   axis = { 0, 1, 0},
   angle = alpha  * 180 / math.pi,
   move = { -hwidth_i, 0, length_i  }
}

box_i = Box{ 
   from = {- hwidth_i, height1,kmin0-100},
   to = { hwidth_i,height2, length_i+ length_it }
}


box_c = Box{ 
   from = {- hwidth_c, height1,length_i + length_it - 1},
   to = { hwidth_c, height2, length_i+ length_it + length_c+1 }
}

l = length_ot / math.cos(beta)
w = l * math.tan(beta)

box_os1 = Transform{ 
   Box{
      from = { 0  , height1, -l },
      to = { w, height2, 0 }
   },
   axis = { 0, 1, 0},
   angle = beta  * 180 / math.pi,
   move = { hwidth_o, 0, length - length_o  }
}

box_os2 = Transform{ 
   Box{
      from = { -w  , height1, -l },
      to = { 0, height2, 0 }
   },
   axis = { 0, 1, 0},
   angle = - beta  * 180 / math.pi,
   move = { -hwidth_o, 0, length - length_o  }
}

box_o = Box{ 
   from = {- hwidth_o, height1, length_i+ length_it + length_c },
   to = { hwidth_o, height2, kmax0+100 }
}

sc_eps:add{ BinaryAndNot{ BinaryAndNot{ box_i, box_is1 }, box_is2 }, depth=1, value=eps_si }
sc_eps:add{ box_c, depth=1, value=eps_si }
sc_eps:add{ BinaryAndNot{ BinaryAndNot{ box_o, box_os1 }, box_os2 }, depth=1, value=eps_si }

--- metal structure!

sc_met = Scene{ value=0 }  


-- .. block

if mcover then

height2 = height2 + height_g

end

box_mcov = Box{ 
   from = { - imax0 - 1, height2, length_i1 },
   to = { imax0 + 1, jmax0 + 1, length - length_o1 }
}

sc_met:add{ box_mcov, depth=1, value=1 }


box_1 = Box{ 
   from = { hwidth_i + width_g, height1, length_i1 },
   to = { imax0 + 1, height2, length_i+ length_it }
}

box_2 = Box{ 
   from = { -hwidth_i - width_g, height1, length_i1 },
   to = { -imax0 - 1, height2, length_i+ length_it }
}

sc_met:add{ box_1, depth=1, value=1 }
sc_met:add{ box_2, depth=1, value=1 }

box_1 = Box{ 
   from = { hwidth_c + width_g, height1, length_i+ length_it },
   to = { imax0 + 1, height2, length_i+ length_it+ length_c }
}

box_2 = Box{ 
   from = { -hwidth_c - width_g, height1, length_i+ length_it },
   to = { -imax0 - 1, height2, length_i+ length_it+ length_c }
}

sc_met:add{ box_1, depth=1, value=1 }
sc_met:add{ box_2, depth=1, value=1 }

box_1 = Box{ 
   from = { hwidth_o + width_g, height1, length - length_o - length_ot },
   to = { imax0 + 1, height2, length - length_o1 }
}

box_2 = Box{ 
   from = { -imax0 - 1, height1, length - length_o - length_ot },
   to = { -hwidth_o - width_g, height2, length - length_o1 }
}

sc_met:add{ box_1, depth=1, value=1 }
sc_met:add{ box_2, depth=1, value=1 }

-- add the missing bits

l = length_it / math.cos(alpha)
w = l * math.tan(alpha)

box_1 = Transform{ 
   Box{
      from = { 0  , height1,  0 },
      to = { w, height2, l }
   },
   axis = { 0, 1, 0},
   angle = - alpha  * 180 / math.pi,
   move = { hwidth_i + width_g, 0, length_i  }
}

box_2 = Transform{ 
   Box{
      from = { -w  , height1,  0 },
      to = { 0 , height2, l }
   },
   axis = { 0, 1, 0},
   angle = alpha  * 180 / math.pi,
   move = { -hwidth_i - width_g, 0, length_i  }
}

sc_met:add{ box_1, depth=1, value=1 }
sc_met:add{ box_2, depth=1, value=1 }

l = length_ot / math.cos(beta)
w = l * math.tan(beta)

box_1 = Transform{ 
   Box{
      from = { 0  , height1, -l },
      to = { w, height2, 0 }
   },
   axis = { 0, 1, 0},
   angle = beta  * 180 / math.pi,
   move = { hwidth_o + width_g, 0, length - length_o  }
}

box_2 = Transform{ 
   Box{
      from = { -w  , height1, -l },
      to = { 0, height2, 0 }
   },
   axis = { 0, 1, 0},
   angle = - beta  * 180 / math.pi,
   move = { -hwidth_o - width_g, 0, length - length_o  }
}

sc_met:add{ box_1, depth=1, value=1 }
sc_met:add{ box_2, depth=1, value=1 }

-- create grids

gj0 = height_bsio2-2
gj1 = height_bsio2+height_bsi+height_wg+2 
gj2 = height_bsio2+height_bsi+height_wg+height_g+2 


grid_eps = Grid{
   from={imin0-1,gj0,kmin0-1},
   to={imax0+1,gj1,kmax0+1} 
}

if mcover then
   gj2 = gj1 +height_g
end

grid_met = Grid{
   from={imin0-1,gj0,kmin0-1},
   to={imax0+1,gj2,kmax0+1} 
}

pad = 80

inj_eps = Grid{
   from = {-hwidth_i-pad,gj0-pad,kil }, 
   to   = {hwidth_i+pad,gj1+pad,kil}
}


prev_eps = Grid{
   yee    = false, 
   from   = {-imax0,gj0,kmin0}, 
   to     = {imax0,gj1,kmax0},
}

prev_met = Grid{
   yee    = false, 
   from   = {-imax0,gj0,kmin0}, 
   to     = {imax0,gj2,kmax0},
}


cfg:CREATE_PREVIEW                  -- create a preview of the scene
{"eps",                             -- filename: preview_"***".vtk
scene=sc_eps,                          -- scene to be added
grid=prev_eps,                     -- grid to be used 
method="default", 
silent=false, 
on=true }

cfg:CREATE_PREVIEW                  -- create a preview of the scene
{"met",                             -- filename: preview_"***".vtk
scene=sc_met,                          -- scene to be added
grid=prev_met,                     -- grid to be used 
method="default", 
silent=false, 
on=true }

cfg:CREATE_GEO                       -- add the scene to the geometry
{"eps",                              -- filename: geo_"***".in 
scene=sc_eps,                           -- scene to be added  
grid=grid_eps,                       -- grid to be used 
method="default",
comps=3,                             -- number of components of permittivity (or permittivity + permeability) 
silent=false,                        -- silient computation, what will happen if we set true?
on=true }                            -- scene will be updated (true)

cfg:CREATE_GEO                       -- add the scene to the geometry
{"met",                              -- filename: geo_"***".in 
scene=sc_met,                        -- scene to be added  
grid=grid_met,                       -- grid to be used 
method="default",
comps=3,                             -- number of components of permittivity (or permittivity + permeability) 
silent=false,                        -- silient computation, what will happen if we set true?
on=true }                            -- scene will be updated (true)


cfg:CREATE_GEO                       -- create the waveguide injection plane
{"inj", 
scene=sc_eps, 
grid=inj_eps, 
method="default",
comps=3, 
silent=false, 
on=true }


--- fire up matlab mode calculator!

cmd = "./luacfg wgmode3k.lua geo_inj.in "..tostring(invwavelength).." "..tostring(betaeff).." 1" 

print(cmd)
os.execute(cmd)


--- GRID Definition

cfg:GRID{

   dim = 3,          
   partition = { 0, 1 },
   ncyc = ncyc,
   dt = dt,  
   irange = { imin0,imax0 },
   jrange = { jmin0,jmax0 },
   krange = { kmin0,kmax0 }

}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{                           
	    { -imax0-1, imax0+1,        1,
	       jmin0, round(height_bsio2), 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_sio2, eps_sio2, eps_sio2 },
	    
            { -imax0-1, imax0+1,        1, 
              round(height_bsio2), jmax0+1, 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_bg, eps_bg, eps_bg }
	 },
	 LOAD_GEO{ "eps" }         
      },
      on = true   
   },

   OUT{
      file = { "VTK", "slice_eps_il" },  
      type = { "Eps", "N" },
      time = { 0, 0, 1 },
      REG{
	 BOX{
	    {  -imax, imax, 1, 
	       jmin, jmax, 1, 
	       kil, kil, 1  }
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_ib" },  
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jc, kib }
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_il" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jc, kl }  
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_ch" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { ich, jc, kc }
	 }
      }
   },

  
   OUT{
      file = { "VTK", "slice_xy_ib" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	       kib, kib, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "slice_xy_il" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, ici, 1, 
	       jmin, jmax, 1, 
	       kil, kil, 1  }
	 }
      }
   },


    OUT{
      file = { "VTK", "slice_xy_i1" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, ici, 1, 
	       jmin, jmax, 1, 
	       ki1, ki1, 1  }
	 }
      }
   },


    OUT{
      file = { "VTK", "slice_xy_i2" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, ici, 1, 
	       jmin, jmax, 1, 
	       ki2, ki2, 1  }
	 }
      }
   },


    OUT{
      file = { "VTK", "slice_xy_c" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, icc, 1, 
	       jmin, jmax, 1, 
	       kc, kc, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "slice_xy_o2" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, ico, 1, 
	       jmin, jmax, 1, 
	       ko2, ko2, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "slice_xy_o1" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  0, ico, 1, 
	       jmin, jmax, 1, 
	       ko1, ko1, 1  }
	 }
      }
   },

     OUT{
       file = { "VTK", "slice_yz_e" },
       type = { "E", "N" },
      time = { 0, ncyc, 250 },
       REG{
	  BOX{
	     { 0, 0, 1, 
	       jmin, jmax, 1,
	       kmin, kmax, 1	       
	    }
	  }
       }
    },

    OUT{
       file = { "VTK", "slice_xz_e" },
       type = { "E", "N" },
      time = { 0, ncyc, 250 },
       REG{
	  BOX{
	     { -imax, imax, 1, 
	       jc, jc, 1,
	       kmin, kmax, 1	       
	    }
	  }
       }
    }

}

--- BOUND Definition (pml boundary condition)

cfg:BOUND{
   config = { bc, 1, 1, 1, 1, 1 },   -- 1 pml 0 means not
   PML{
      cells = 11,                   -- no of pml cells
      pot = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- MAT Definition(s) 

cfg:MAT{
   PEC{},
--   DRUDE{
--      invlambdapl = 0.0471404520791,
--      gammapl = 1.e-3,
--      order = 2
--   },
   REG{
      LOAD_GEO{ "met" }
   }
}

--- SRC Definition(s) (source)

cfg:SRC{
   TFSFINJ{                             -- total field scattered field injection 
      invlambda = invwavelength,        -- inverse of wavelength (frequency)
      amplitude = 1.,                   -- amplitude of the source
      pulse = { 
	 shape   = "Gaussian",          -- gaussian envelope
	 width   = pulsehwhm,           -- pulsehwhm
	 offset  = 0,                    
         attack  = pulsehsteps,        
         sustain = 0, 
         decay   = pulsehsteps   
      },
      planewave = injplane  -- this gives the position of injection plane and refractive index of injection plane
   },
   REG{
      LOAD{ injfile }   -- input waveguide mode
   },
   on = true
}


--- PSPEC


cfg:DIAG{
   PSPEC{
      file = "sib",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ici, isi, jc1, jc3, js, kib, kib, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "si1",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ici, isi, jc1, jc3, js, ki1, ki1, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "si2",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ici, isi, jc1, jc3, js, ki2, ki2, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "sc",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, icc, isc, jc1, jc3, js, kc, kc, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "so2",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ico, iso, jc1, jc3, js, ko2, ko2, 1 } 
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "so1",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ico, iso, jc1, jc3, js, ko1, ko1, 1 } 
      }
   }
}



cfg:CREATE()