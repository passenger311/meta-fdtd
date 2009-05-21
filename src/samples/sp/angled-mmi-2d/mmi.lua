cfg = CONFIG{scenes=true}  -- scences takes either true or false 

TASKID = select(1,...)

dofile("scale.lua")  -- geometrical and material parameters are defined in scale file

-- grid dimensions: i,j,k

imax = round( hwidth )
jmax = round( hlength )  -- propagation direction
kmax = round( height )

--- (print("computational window without pml: ", imin, imax, jmin, jmax, kmin, kmax))

imin = -imax
jmin = - jmax
kmin = 0

print("irange (grid) = ", imin, imax)
print("jrange (grid) = ", jmin, jmax)
print("krange (grid) = ", kmin, kmax)
print("ncyc          = ", ncyc)

print("irange (real) = ", imin*real_dx, imax*real_dx)
print("jrange (real) = ", jmin*real_dx, jmax*real_dx)
print("krange (real) = ", kmin*real_dx, kmax*real_dx)


--- full domain including PMLs

pad1 = 0

imin0 = imin -pad1 - cpml
imax0 = imax +pad1 + cpml
jmin0 = jmin -pad1 - cpml
jmax0 = jmax +pad1 + cpml
kmin0 = kmin -pad1 - cpml
kmax0 = kmax +pad1 + cpml



--- create dielectric structure (scenes)

mmi = Scene{ value=eps_bg }   -- value represents dielectric constant of the background value=1 means air background

---wg = Scene{value=1.}
--- cladding block it does not matter if the block is larger than its actual size (-/+ 100)
box_bsio2 = Box{ 
   from={-imax0-100,jmin0-100,kmin0-100},
   to={imax0+100,jmax0+100,height_bsio2}
}

--- mmi block (silicon)

box_bsi = Box{ 
   from={-imax0-100,jmin0-100,height_bsio2},   --  height_bsio2 represents jmin0 for this block
   to={imax0+100,jmax0+100,height_bsio2+height_bsi} -- height_bsio2+height_bsi gives height of silicon block
}

mmi:add{ box_bsio2, depth=1, value=eps_sio2 }  -- paint geometrical objects on the scene depth represents the permittivity
mmi:add{ box_bsi, depth=1, value=eps_si }      -- lower value of depth of the object is painted on the top

mmwg =  Transform{ 
   Box{ 
      from = { -hwidth_mmi, -hlength_mmi, ys },
      to = { hwidth_mmi, hlength_mmi, ye }
   },
   axis = { 0, 0, 1 },
   angle = -angle
}

inch = Box{ 
   from = { x8, z8-100, ys },
   to = { x7, z7+5, ye }
}

mmi:add{ mmwg, depth = 1, value=eps_si }
mmi:add{ inch, depth = 1, value=eps_si }

from = { x10, z9-5, ys }
to = { x9, z10+100, ye }

for i = 1, numchannel do

   outch = Box{ 
      from = from,
      to = to
   }

   mmi:add{ outch, depth = 1, value=eps_si }

   from = { from[1]-ch_dx, from[2]-ch_dz, from[3] }
   to = { to[1]-ch_dx, to[2], to[3] }

end

gkmin = height_bsio2                         -- cladding height          
gkmax = height_bsio2+height_bsi+height_wg    -- total height of the structure

--- specify grid for the scene, only real structure will be discretised

grid_eps = Grid{
   from={imin0,jmin0,kc-1},  -- left corner of the geometry pml not included 
   to={imax0,jmax0,kc+1}   -- upper front right corner of the geometry pml not included
}

pad = 80
-- waveguide grid in injection plane
grid_inj = Grid{
   from = {x8-pad, jinj, gkmin-pad }, 
   to   = {x8+2*hwidth_wg+pad, jinj ,gkmax+pad}
}

print("CENTER OF INJ I/J: ", x8+hwidth_wg, jinj)

-- coarse grid for .VTK-preview
grid_prev = Grid{
   yee    = false,                            -- single permittivity component only
   from={imin0,jmin0,kc-1},  -- left corner of the geometry pml not included 
   to={imax0,jmax0,kc+1}   -- upper front right corner of the geometry pml not included
}


cfg:CREATE_PREVIEW   -- create a preview of the scene
{"mmi",              -- filename: preview_"***".vtk
scene=mmi,           -- scene to be added
grid=grid_prev,      -- grid to be used 
method="default", 
silent=false, 
on=true }

cfg:CREATE_GEO   -- add the scene to the geometry
{"mmi",          -- filename: geo_"***".in 
scene=mmi,       -- scene to be added  
grid=grid_eps,   -- grid to be used 
method="default",
comps=3,         -- number of components of permittivity (or permittivity + permeability) 
silent=false,    -- silient computation, what will happen if we set true?
on=true }        -- scene will be updated (true)

cfg:CREATE_GEO   -- create the waveguide injection plane
{"inj", 
scene=mmi, 
grid=grid_inj, 
method="default",
comps=3, 
silent=false, 
on=true }

--- fire up matlab mode calculator!

print("> luacfg mode.lua "..tostring(invwavelength).." "..tostring(betaeff).." 1")
os.execute("luacfg mode.lua "..tostring(invwavelength).." "..tostring(betaeff).." 1")


--- GRID Definition

cfg:GRID{

   dim = 2,                  -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncyc,              -- total number simulations time (real total time = ncyc*dt)
   dt = dt,                  -- time step 
   irange = { imin0,imax0 }, -- range of computational cells in x direction 
   jrange = { jmin0,jmax0 }, -- range of computational cells in y direction
   krange = { kc,kc }  -- range of computational cells in z direction
}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{                            -- background block
	    { 
	       imin0, imax0+1,        1,  -- coordinate start, end and step 
	       jmin0, jmax0+1,        1,
	       kc, kc,                1, ":", 
	       eps_bg, eps_bg, eps_bg 
	    }
	 },
	 LOAD_GEO{ "mmi" }          -- will load mmi geometry 
      },
      on = true   
   },

   OUT{
      file = { "GPL", "point_e_inb" },  -- gpl (gnuplot format) electric field before source
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },              -- start end step 
      REG{
	 POINT{ 
	    { 0, jinb }     -- this is the point where we calculate the field evloution with time
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_inf" },   -- just after source 
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jinf }  
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_mid" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, 0, kc }   -- middle point of the structure  
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_ch1" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { ich[1], jch[1] }   -- end point (at the output channel)
	 }
      }
   },

   OUT{
      file = { "VTK", "slice_e_jinb" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jinb, jinb, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "slice_e_jinf" },
      type = { "E", "N" },
      time = { 0, ncyc, 250 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jinf, jinf, 1  }
	 }
      }
   },

    OUT{
       file = { "VTK", "slice_e_xz" },
       type = { "E", "N" },
      time = { 0, ncyc, 250 },
       REG{
	  BOX{
	     { imin, imax, 1, 
	       jmin, jmax, 1	       
	    }
	  }
       }
    }

}

--- BOUND Definition (pml boundary condition)

cfg:BOUND{
   config = { 1, 1, 1, 1, 1, 1 },   -- 1 pml 0 means not
   PML{
      cells = 11,                   -- no of pml cells
      pot = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
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
      planewave = injplane
   },
   REG{
      LOAD{ injfile }                   -- input waveguide mode
   },
   on = true
}


--- PSPEC

cfg:DIAG{
   PSPEC{
      file = "inb1",
      time = {1, 4*pulsehsteps, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { iin-idelta, iin+idelta, istep, jinb, jinb, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "inf1",
      time = {1, 4*pulsehsteps, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { iin-idelta, iin+idelta, istep, jinf, jinf, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "inb",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { iin-idelta, iin+idelta, istep, jinb, jinb, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "inf",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { iin-idelta, iin+idelta, istep, jinf, jinf, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "in0",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { iin-idelta, iin+idelta, istep, jin0, jin0, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "ch1",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { ich[1]-idelta, ich[1]+idelta, istep, jch[1], jch[1], 1 }
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "ch2",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { ich[2]-idelta, ich[2]+idelta, istep, jch[2], jch[2], 1 }
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "ch3",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
	 { ich[3]-idelta, ich[3]+idelta, istep, jch[3], jch[3], 1 }
      }
   }
}



--- CREATE: config.<part>.in
cfg:CREATE()