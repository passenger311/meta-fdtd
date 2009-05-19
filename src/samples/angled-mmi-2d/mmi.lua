cfg = CONFIG{scenes=true}  -- scences takes either true or false 

taskid = select(1,...)

dofile("scale.lua")  -- geometrical and material parameters are defined in scale file

-- grid dimensions: i,j,k

imax = round( hwidth )
jmax = round( height )         -- j along z axis 
kmax = round( hlength )        -- k along y axis this is the direction where wave propagates

--- (print("computational window without pml: ", imin, imax, jmin, jmax, kmin, kmax))

imin = -imax
jmin = 0
kmin = - kmax

print("irange (grid) = ", imin, imax)
print("jrange (grid) = ", jmin, jmax)
print("krange (grid) = ", kmin, kmax)

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
   to={imax0+100,height_bsio2,kmax0+100}
}

--- mmi block (silicon) 
box_bsi = Box{ 
   from={-imax0-100,height_bsio2,kmin0-100},   --  height_bsio2 represents jmin0 for this block
   to={imax0+100,height_bsio2+height_bsi,kmax0+100} -- height_bsio2+height_bsi gives height of silicon block
}

mmi:add{ box_bsio2, depth=1, value=eps_sio2 }  -- paint geometrical objects on the scene depth represents the permittivity
mmi:add{ box_bsi, depth=1, value=eps_si }      -- lower value of depth of the object is painted on the top


mmwg =  Transform{ 
   Box{ 
      from = { -hwidth_mmi, ys, -hlength_mmi },
      to = { hwidth_mmi, ye, hlength_mmi }
   },
   axis = { 0, 1, 0 },
   angle = angle
}

inch = Box{ 
   from = { x8, ys, z8-100 },
   to = { x7, ye, z7+5 }
}

mmi:add{ mmwg, depth = 1, value=eps_si }
mmi:add{ inch, depth = 1, value=eps_si }

from = { x10, ys, z9-5 }
to = { x9, ye, z10+100 }

for i = 1, numchannel do

   outch = Box{ 
      from = from,
      to = to
   }

   mmi:add{ outch, depth = 1, value=eps_si }

   from = { from[1]-ch_dx, from[2], from[3]-ch_dz }
   to = { to[1]-ch_dx, to[2], to[3] }

end

gjmin = height_bsio2                           -- cladding height          
gjmax = height_bsio2+height_bsi+height_wg      -- total height of the structure

--- specify grid for the scene, only real structure will be discretised

grid_eps = Grid{
   from={imin0,gjmin-10,kmin0},  -- left corner of the geometry pml not included 
   to={imax0,gjmax+10,kmax0}   -- upper front right corner of the geometry pml not included
}

pad = 80
-- waveguide grid in injection plane
grid_inj = Grid{
   from = {x8-pad,gjmin-pad,kinj }, 
   to   = {x8+2*hwidth_wg+pad,gjmax+pad,kinj}
}
-- coarse grid for .VTK-preview
grid_prev = Grid{
   yee    = false,                            -- single permittivity component only
   from={imin0,gjmin-10,kmin0},  -- left corner of the geometry pml not included 
   to={imax0,gjmax+10,kmax0}   -- upper front right corner of the geometry pml not included
--   from   = {imin0,gjmin-10,kmin0},          -- lower rear left corner of scene-object 
--   to     = {imax0,gjmax+10,kmax0},           -- upper front right corner of scene-object
--   offset = {imin0,jmin0,kmin0 },             -- this is origin of coordinate system, not really offset
--   cells  = {gjmax+10-gjmin+10, 2*hwidth , 2*hlength}  -- this gives number of cells in each direction to grid scene
}

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

cfg:CREATE_PREVIEW   -- create a preview of the scene
{"mmi",              -- filename: preview_"***".vtk
scene=mmi,           -- scene to be added
grid=grid_prev,      -- grid to be used 
method="default", 
silent=false, 
on=true }

--- GRID Definition

cfg:GRID{

   dim = 3,                  -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncyc,              -- total number simulations time (real total time = ncyc*dt)
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
	    { 
	       imin0, imax0+1,        1,  -- coordinate start, end and step 
	       jmin0, round(height_bsio2)-1, 1, 
	       kmin0, kmax0+1,        1, ":", 
	       eps_sio2, eps_sio2, eps_sio2 
	    },   -- dielectric constant  
	    
            { 
	       imin0, imax0+1,        1, 
	       round(height_bsio2), jmax0+1, 1, 
	       kmin0, kmax0+1,        1, ":", 
	       eps_bg, eps_bg, eps_bg 
	    }
	 },
	 LOAD_GEO{ "mmi" }          -- will load mmi geometry 
      },
      on = true   
   },


   OUT{
      file = { "SET", "mmi_eps" },   -- set here represents file format (mmi_eps.set)
      type = { "Eps", "N" },
      time = { 0, 0, 1 },            -- start end and step
      REG{
	 BOX{
	    {  imin,imax, 1, 
	       jmin,jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_injb" },  -- gpl (gnuplot format) electric field before source
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },              -- start end step 
      REG{
	 POINT{ 
	    { 0, jc, kinj-1 }     -- this is the point where we calculate the field evloution with time
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_injf" },   -- just after source 
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jc, kinj+1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_mid" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, jc, 0 }   -- middle point of the structure  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_end" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { iout1, jc, kfft3 }   -- end point (at the output channel)
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft1" },   -- energy density at inch (before source!)
      type = { "En", "S", ".F." },             -- S stands for sum En stands for six components of E and H fields
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
           { iin-deltai, iin+deltai, 2, jc, jc, 1, kfft1, kfft1, 1 }  -- this defines range of points where we collect En and sum it   
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft2" },   -- energy density  at inch (after source!)
      type = { "En", "S", ".F." },             -- S stands for sum En stands for six components of E and H fields
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
           { iin-deltai, iin+deltai, 2, jc, jc, 1, kfft2, kfft2, 1 }  -- this defines range of points where we collect En and sum it   
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_out1" },   -- energy density at outch1
      type = { "En", "S", ".F." },             -- S stands for sum En stands for six components of E and H fields
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
           { iout1-deltai, iout1+deltai, 2, jc, jc, 1, kfft3, kfft3, 1 }  -- this defines range of points where we collect En and sum it   
	 }
      }
   },

  
   OUT{
      file = { "VTK", "mmi_slice1_xy_eps" },
      type = { "Eps", "N" },
      time = { 0, 0, 1 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "mmi_slice0_xy_e" },
      type = { "E", "N" },
      time = { 0, ncyc, 200 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	       kinj-1, kinj-1, 1  }
	 }
      }
   },
   
    OUT{
       file = { "VTK", "mmi_slice1_xz_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 100 },
       REG{
	  BOX{
	     { imin, imax, 1, 
	       jc, jc, 1,
	       kmin, kmax, 1	       
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
      -- set psi 0./90. for TE/ TM
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=neff }  -- this gives the position of injection plane and refractive index of injection plane
   },
   REG{
      LOAD{ "tfsf90.set" }   -- input waveguide mode
   },
   on = true
}


--- PSPEC


cfg:DIAG{
   PSPEC{
      file = "inb",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { iin-deltai, iin+deltai, 2, jc1, jc2, 2, kfft1, kfft1, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "inf",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { iin-deltai, iin+deltai, 2, jc1, jc2, 2, kfft2, kfft2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "out1",
      time = {1, ncyc, 1},
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { iout1-deltai, iout1+deltai, 2, jc, jc, 2, kfft3, kfft3, 1 }
      }
   }
}



--- CREATE: config.<part>.in
cfg:CREATE()