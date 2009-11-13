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
kmin0 = kmin -pad1          -- cpml is included in the structure size
kmax0 = kmax +pad1       -- cpml is included in the structure size


--- dielectric structure!

sc_eps = Scene{ value=eps_bg }  

-- ... sio2 base

box_bsio2 = Box{
   from={-imax0-100,jmin0-100,kmin0-100},
   to  ={imax0+100,height_bsio2,kmax0+100}
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

if cladding then

--- cladding

-- .. block

box_1 = Box{ 
   from = { hwidth_i + width_gi, height1, length_i1 },
   to = { hwidth_i + width_gi + width_ii, height2, length_i }
}

box_2 = Box{ 
   from = { -hwidth_i - width_gi, height1, length_i1 },
   to = { -hwidth_i - width_gi - width_ii, height2, length_i }
}

sc_eps:add{ box_1, depth=1, value=eps_si }
sc_eps:add{ box_2, depth=1, value=eps_si }

box_1 = Box{ 
   from = { hwidth_c + width_gc, height1, length_i+ length_it },
   to = { hwidth_c + width_gc + width_cc, height2, length_i+ length_it+ length_c }
}

box_2 = Box{ 
   from = { -hwidth_c - width_gc, height1, length_i+ length_it },
   to = { -hwidth_c - width_gc - width_cc, height2, length_i+ length_it+ length_c }
}

sc_eps:add{ box_1, depth=1, value=eps_si }
sc_eps:add{ box_2, depth=1, value=eps_si }

box_1 = Box{ 
   from = { hwidth_o + width_go, height1, length - length_o },
   to = { hwidth_o + width_go + width_oo, height2, length - length_o1 }
}

box_2 = Box{ 
   from = { -hwidth_o - width_go, height1, length - length_o },
   to = { -hwidth_o - width_go - width_oo, height2, length - length_o1 }
}

sc_eps:add{ box_1, depth=1, value=eps_si }
sc_eps:add{ box_2, depth=1, value=eps_si }


-- add the missing bits

l = length_it / math.cos(alpha)
w = l * math.tan(alpha)

xs = math.min(hwidth_i + width_gi, hwidth_c + width_gc)

box_1 = Box{ from = { xs, height2, length_i },
	     to = { infty, height1, length_i+ length_it } }


box_2 = Box{ from = { -xs, height2, length_i },
	     to = { -infty, height1, length_i+ length_it } }


box_1a =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { -infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = - alpha  * 180 / math.pi,
   move = { hwidth_i + width_gi , 0, length_i  }
}

box_1b =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = - alpha2  * 180 / math.pi,
   move = { hwidth_i + width_gi + width_ii, 0, length_i  }
}


box_2a =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = alpha  * 180 / math.pi,
   move = { -hwidth_i - width_gi , 0, length_i  }
}

box_2b =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { -infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle =  alpha2  * 180 / math.pi,
   move = { -hwidth_i - width_gi - width_ii, 0, length_i  }
}


sc_eps:add{ BinaryAndNot{BinaryAndNot{box_1,box_1b},box_1a}, depth=1, value=eps_si }
sc_eps:add{ BinaryAndNot{BinaryAndNot{box_2,box_2b},box_2a}, depth=1, value=eps_si }

xs = math.min(hwidth_c + width_gc, hwidth_o + width_go)

box_1 = Box{ from = { xs, height2, length_i+ length_it+ length_c },
	     to = { infty, height1, length - length_o } }

box_2 = Box{ from = { -xs, height2, length_i+ length_it+ length_c },
	     to = { -infty , height1, length - length_o } }



box_1a =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { -infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = beta  * 180 / math.pi,
   move = { hwidth_c + width_gc , 0, length - length_o - length_ot  }
}

box_1b =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = beta2  * 180 / math.pi,
   move = { hwidth_c + width_gc + width_cc, 0, length - length_o - length_ot  }
}

box_2a =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = - beta  * 180 / math.pi,
   move = { -hwidth_c - width_gc , 0, length - length_o - length_ot  }
}

box_2b =  Transform{ 
   Box{
      from = { 0  , height1,  -infty },
      to = { - infty, height2, infty }
   },
   axis = { 0, 1, 0},
   angle = - beta2  * 180 / math.pi,
   move = {  -hwidth_c - width_gc - width_cc, 0, length - length_o - length_ot  }
}

sc_eps:add{ BinaryAndNot{BinaryAndNot{box_1,box_1b},box_1a}, depth=1, value=eps_si }
sc_eps:add{ BinaryAndNot{BinaryAndNot{box_2,box_2b},box_2a}, depth=1, value=eps_si }

end

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


cfg:CREATE_PREVIEW                  -- create a preview of the scene
{"eps",                             -- filename: preview_"***".vtk
scene=sc_eps,                          -- scene to be added
grid=prev_eps,                     -- grid to be used 
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
	    { ib, jc, kib }            -- 0 was there before in the position of ib
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_il" },  -- record time response of field
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { ib, jc, kil }    -- k i and letter L
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

  --[[
  OUT{
     file = { "VTK", "slice_xy_ib" },
     type = { "E", "N" },
     time = { 1000, ncyc, 1000 },
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
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    {  0, ici, 1, 
	       jmin, jmax, 1, 
	       kil, kil, 1  }
	 }
      }
   }, --]]


OUT{
      file = { "VTK", "slice_xy_i11" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    --[[{  0, iinj, 1, 
	       jinj1, jinj2, 1, 
	       kinj, kinj, 1  }--]]
              
               {-iinj-pad, iinj+pad, 1, 
                 jinj1-1, jinj2+1, 1, 
                 kil, kil, 1} 
	 }
      }
   },

--[[
OUT{
      file = { "VTK", "slice_xy_i2" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
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
      time = { 1000, ncyc, 1000 },
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
      time = { 1000, ncyc, 1000 },
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
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    {  0, ico, 1, 
	       jmin, jmax, 1, 
	       ko1, ko1, 1  }
	 }
      }
   },
--]]
     
OUT{
       file = { "VTK", "slice_yz_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
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
      time = { 1000, ncyc, 1000 },
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
      pot   = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- MAT Definition(s) 

--cfg:MAT{
--   PEC{},
--   DRUDE{
--      invlambdapl = 0.0471404520791,
--      gammapl = 1.e-3,
--      order = 2
--   },
--   REG{
--      LOAD_GEO{ "met" }
--   }
--}



--- SRC Definition(s) (source)


cfg:SRC{
   TFSFINJ{                             -- total field scattered field injection 
      invlambda = invwavelength,        -- inverse of wavelength (frequency)
      amplitude = 1.,                   -- amplitude of the source
      pulse = pulse,
      planewave = { phi=0, theta=0, psi=90, nrefr=betaeff}
   },
   REG{
--      LOAD{"tfsfex.set"}
--      BOX{ { -iinj-1, iinj+1, 1, jinj1-1, jinj2+1, 1, kinj, kinj, 1} } 
        BOX{ 
	{ -iinj, iinj, 1, jinj1-1, jinj2+1, 1, kil, kil, 1} 
	}
   },
   on = true
}

--- PSPEC


cfg:DIAG{
   PSPEC{
      file = "sib",
      time = {1, ncyc, 10},
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
      time = {1, ncyc, 10},
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
      time = {1, ncyc, 10},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, ici, isi, jc1, jc3, js, ki2, ki2, 1 } 
      }
   }
}

---same width as output so2
cfg:DIAG{
   PSPEC{
      file = "si22",
      time = {1, ncyc, 10},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
         {0, ico, iso, jc1, jc3, js, ko2, ko2, 1}
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "sc",
      time = {1, ncyc, 10},
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
      time = {1, ncyc, 10},
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
      time = {1, ncyc, 10},
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