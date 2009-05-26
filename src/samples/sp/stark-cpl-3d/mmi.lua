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


--- create dielectric structure (scenes)

mmi = Scene{ value=eps_bg }  

-- do sio2 base

box_bsio2 = Box{
   from={-imax0-100,jmin0-100,kmin0-100},
   to={imax0+100,height_bsio2,kmax0+100}
}

mmi:add{ box_bsio2, depth=1, value=eps_sio2 }

-- do silicon base

box_bsi = Box{ 
   from={-imax0-100,height_bsio2,kmin0-100},
   to={imax0+100,height_bsio2+height_bsi,kmax0+100}
}

mmi:add{ box_bsi, depth=1, value=eps_si }     

joffs1 = height_bsio2+height_bsi-1
joffs2 = height_bsio2 + height_bsi + height_wg

box_wg1 = Box{
   from = { -hwidth_wg, joffs1, kmin0-100 },
   to = { hwidth_wg, joffs2, length_wg1+1 }
}

mmi:add{ box_wg1, depth=1, value=eps_si }

box_mmi = Box{
   from = { -hwidth_mmi, joffs1, length_wg1 },
   to = { hwidth_mmi, joffs2, length_wg1+length_mmi }
}

mmi:add{ box_mmi, depth=1, value=eps_si }

box_wg21 = Transform{  
   Box{
      from = { 0, joffs1, -2 },
      to = { 2*hwidth_wg, joffs2, length_wg2+100000 }
   },
   axis = { 0, 1, 0 },
   angle = angle,
   move = { p8[1],0,p8[2] }  
}

box_wg22 = Transform{  
   Box{
      from = { -2*hwidth_wg, joffs1, -2 },
      to = { 0, joffs2, length_wg2+100000 }
   },
   axis = { 0, 1, 0 },
   angle = -angle,
   move = { -p8[1],0,p8[2] }  
}


box_wg21_corr = BinaryAndNot{
   box_wg21,
   Box{
      from = { imin0-100, joffs1, length_wg1-100 },
      to = { hwidth_mmi, joffs2, length_wg1 }
   }
  
}
box_wg22_corr = BinaryAndNot{
   box_wg22,
   Box{
      from = { imin0-100, joffs1, length_wg1-100 },
      to = { hwidth_mmi, joffs2, length_wg1 }
   }
  
}

mmi:add{ box_wg21_corr, depth=1, value=eps_si }
mmi:add{ box_wg22_corr, depth=1, value=eps_si }


-- the vertical layer with the structure embedded

gjmin = height_bsio2
gjmax = height_bsio2+height_bsi+height_wg 


grid_eps = Grid{
   from={-imax0-1,gjmin-10,kmin0-1},
   to={imax0+1,gjmax+10,kmax0+1} 
}

pad = 80

grid_inj = Grid{
   from = {-hwidth_wg-pad,gjmin-pad,kinj }, 
   to   = {hwidth_wg+pad,gjmax+pad,kinj}
}

grid_prev = Grid{
   yee    = false, 
   from   = {-imax0,gjmin-10,kmin0}, 
   to     = {imax0,gjmax+10,kmax0},
}


cfg:CREATE_PREVIEW                  -- create a preview of the scene
{"mmi",                             -- filename: preview_"***".vtk
scene=mmi,                          -- scene to be added
grid=grid_prev,                     -- grid to be used 
method="default", 
silent=false, 
on=true }

cfg:CREATE_GEO                       -- add the scene to the geometry
{"mmi",                              -- filename: geo_"***".in 
scene=mmi,                           -- scene to be added  
grid=grid_eps,                       -- grid to be used 
method="default",
comps=3,                             -- number of components of permittivity (or permittivity + permeability) 
silent=false,                        -- silient computation, what will happen if we set true?
on=true }                            -- scene will be updated (true)

cfg:CREATE_GEO                       -- create the waveguide injection plane
{"inj", 
scene=mmi, 
grid=grid_inj, 
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

   dim = 3,                  -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncyc,              -- simulation time (real total time = ncyc*dt)
   dt = dt,                  -- time step 
   irange = { imin0,imax0 }, -- range of computational cells in x direction 
   jrange = { jmin0,jmax0 }, -- range of computational cells in y direction
   krange = { kmin0,kmax0 }  -- range of computational cells in z direction

}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{                            -- background block
	    { -imax0-1, imax0+1,        1,  -- coordinate start, end and step 
	       jmin0, round(height_bsio2), 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_sio2, eps_sio2, eps_sio2 },   -- dielectric constant  
	    
            { -imax0-1, imax0+1,        1, 
              round(height_bsio2), jmax0+1, 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_bg, eps_bg, eps_bg }
	 },
	 LOAD_GEO{ "mmi" }          -- will load mmi geometry 
      },
      on = true   
   },

   OUT{
      file = { "VTK", "slice_eps_inj" },  
      type = { "Eps", "N" },
      time = { 0, 0, 1 },                 -- start end step 
      REG{
	 BOX{
	    {  -imax, imax, 1, 
	       jmin, jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_injb" },  
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },                 -- start end step 
      REG{
	 POINT{ 
	    { 0, jc, kinb }                -- this is the point where we calculate the field evloution with time
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_injf" },   -- after source 
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jc, kinf }  
	 }
      }
   },

   OUT{
      file = { "GPL", "point_e_ch" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { ich, jc, kch }   -- end point (at the output channel)
	 }
      }
   },

  
   OUT{
      file = { "VTK", "slice0_xy_inb" },
      type = { "E", "N" },
      time = { 0, 1000, 250 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	       kinb, kinb, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "slice0_xy_inj" },
      type = { "E", "N" },
      time = { 0, 1000, 250 },
      REG{
	 BOX{
	    {  -imax, imax, 1, 
	       jmin, jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },
   
    OUT{
       file = { "VTK", "slice1_xz_e" },
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
      file = "inb",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = injplane
   },
   REG{
      BOX{ 
        { 0, idelta, istep, jc1, jc2, jstep, kinb, kinb, 1 } 
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
        { 0, idelta, istep, jc1, jc2, jstep, kinf, kinf, 1 } 
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
        { 0, idelta, istep, jc1, jc2, jstep, kin0, kin0, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "ch1",
        time = { 1, ncyc, 1 },
        mode = "S",
        polarize = { phi=injplane.phi, theta=angle, 
		     psi=injplane.psi, nrefr=injplane.nrefr }
   },
   REG{
      BOX{ 
        { 0, idelta, istep, jc1, jc2, jstep, kin0, kin0, 1 } 
      }
   }
}
--- CREATE: config.<part>.in
cfg:CREATE()