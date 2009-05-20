cfg = CONFIG{scenes=true}  -- scences takes either true or false 

taskid = select(1,...)

dofile("scale.lua")  -- geometrical and material parameters are defined in scale file

-- grid dimensions: i,j,k

width = math.floor((2*real_width_wg+real_width_mmi)/real_dx) -- width in grid points width of mmi plus twice wg width

imax = math.floor( (width - 1) / 2  )   -- i along x axis
jmax = math.floor( height - 1 )         -- j along z axis 
kmax = math.floor( length - 1 )         -- k along y axis this is the direction where wave propagates

--- (print("computational window without pml: ", imin, imax, jmin, jmax, kmin, kmax))

print("irange (grid) = ", -imax, imax)
print("jrange (grid) = ", 0,     jmax)
print("krange (grid) = ", 0,     kmax)

print("irange (real) = ", -imax*real_dx, imax*real_dx)
print("jrange (real) = ", 0,             jmax*real_dx)
print("krange (real) = ", 0,             kmax*real_dx)



--- fix some parameters

imin = 0
jmin = 0
kmin = 0

imin0 = imin 
imax0 = imax + cpml
jmin0 = jmin - cpml
jmax0 = jmax + cpml
kmin0 = kmin - cpml
kmax0 = kmax + cpml


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

dofile("mmi_structure.lua")                    -- this file mmi_structure where we defined rotation of channels 

gjmin = height_bsio2                           -- cladding height          
gjmax = height_bsio2+height_bsi+height_wg      -- total height of the structure

--- specify grid for the scene, only real structure will be discretised

grid_eps = Grid{
   from={-imax,gjmin-10,0},  -- left corner of the geometry pml not included 
   to={imax,gjmax+10,kmax}   -- upper front right corner of the geometry pml not included
}

pad = 80
-- waveguide grid in injection plane
grid_inj = Grid{
   from={-hwidth_wg-pad,gjmin-pad,kinj }, 
   to={hwidth_wg+pad,gjmax+pad,kinj}
}
-- coarse grid for .VTK-preview
grid_prev =  Grid{
   yee=false,                    -- single permittivity component only
   from={-imax0,gjmin-10,kmin0}, -- lower rear left corner of scene-object 
   to={imax0,gjmax+10,kmax0},    -- upper front right corner of scene-object
   offset={-imax0,0,0},          -- this is origin not really offset
   cells={70,150,150}            -- this gives number of cells in each direction to grid scene
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
	    { imin0, imax0+1,        1,  -- coordinate start, end and step 
              jmin0, height_bsio2-1, 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_sio2, eps_sio2, eps_sio2 },   -- dielectric constant  
	    
            { imin0, imax0+1,        1, 
              height_bsio2, jmax0+1, 1, 
              kmin0, kmax0+1,        1, ":", 
              eps_bg, eps_bg, eps_bg }
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
	    { 0, yc, kinj-1 }     -- this is the point where we calculate the field evloution with time
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_injf" },   -- after source 
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, yc, kinj+1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_mid" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 1 },
      REG{
	 POINT{ 
	    { 0, yc, math.floor(kmax/2) }   -- middle point of the structure  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_end" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { hwidth_sep, yc, kmax }   -- end point (at the output channel)
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft1" },   -- energy density  
      type = { "En", "S", ".F." },             -- S stands for sum En stands for six components of E and H fields
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	 --   { 1, 21, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 }  
           { imin0, hwidth_mmi, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 }  -- this defines range of points where we collect En and sum it   
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft2" },   
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    --{ hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 } 
           { imin0, hwidth_mmi, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }     
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft3" }, -- energy density
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    -- { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft3, kfft3, 1 }  
           { imin0, hwidth_mmi, 3, yc-10, yc+10, 3, kfft3, kfft3, 1 }    
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
      file = { "VTK", "mmi_slice1_xy_e" },
      type = { "E", "N" },
      time = { 0, ncyc, 200  },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "mmi_slice2_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 200 },
      REG{
	 BOX{
	    {  imin, imax,  1, 
	       jmin, jmax,  1, 
	      kfft1, kfft1, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "mmi_slice3_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 200 },
      REG{
	 BOX{
	    { imin, imax, 1, 
	      jmin, jmax, 1, 
	      math.floor(length_wg1/2),math.floor(length_wg1/2) , 1  }
	 }
      }
   },

    OUT{
       file = { "VTK", "mmi_slice4_xy_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 200 },
       REG{
	  BOX{
	     { imin, imax, 1, 
	       jmin, jmax, 1, 
	       kfft2,kfft2, 1  }
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
	       yc, yc, 1,
	       kmin, kmax, 1	       
	    }
	  }
       }
    }

}

--- BOUND Definition (pml boundary condition)

cfg:BOUND{
   config = { 0, 1, 1, 1, 1, 1 },   -- 1 pml 0 means not
   PML{
      cells = 11,                   -- no of pml cells
      pot = 3.2,                     
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s) (source)

cfg:SRC{
   TFSFINJ{                         -- total field scattered field injection 
      invlambda = invwavelength,    -- inverse of wavelength (frequency)
      amplitude = 1.,               -- amplitude of the source
      pulse = { 
	 shape="Gaussian",          -- gaussian envelope
	 width=pulsehwhm,           -- pulsehwhm
	 offset=0,                    
         attack=pulsehsteps,        
         sustain=0, 
         decay=pulsehsteps   
      },
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=neff }  -- this gives the position of injection plane and refractive index of injection plane
   },
   REG{
      LOAD{ "tfsfinj.set" }   -- input waveguide mode
   },
   on = true
}


--- PSPEC


cfg:DIAG{
   PSPEC{
      file = "mmi1",
      time = { 1, ncyc-128, 1 },  
  --  time = {1, ncyc, 1},
    mode = "Eap",
 --   mode = "S",
      phasewrap = { 1, 0 },
   -- polarize = { phi=0, theta=0, psi=90.0 }
      polarize = { phi=0, theta=30.0, psi=0.0 }
   },
   REG{
      BOX{ 
--	 { 0, imax-1, 3, yc, yc, 1, kfft1, kfft1, 1 }
        { imin, imax-1, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 } 
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "mmi2",
  --    time = { 128, ncyc, 1 },
    time = { 1, ncyc, 1 },
 --     reffile = "mmi1",
      phasewrap = { 1, 0 },
      mode = "Eap",
--    mode = "S",
 --   polarize = { phi=0, theta=0, psi=90.0 }
      polarize = { phi=0, theta=30.0, psi=0.0 }
   },
   REG{
      BOX{ 
	 --{ 0, imax-1, 3, yc, yc, 1, kfft3, kfft3, 1 } 
           { imin, imax-1, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }   
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "mmi3",
  --    time = { 128, ncyc, 1 },
      time = { 1, ncyc, 1 },
  --    reffile = "mmi1",
      phasewrap = { 1, 0 },
   --    mode = "S",
      mode = "Eap",
  --    polarize = { phi=0, theta=0, psi=90.0 }
      polarize = { phi=0, theta=30.0, psi=0.0 }
   },
   REG{
      BOX{ 
	 --{ 0, imax-1, 3, yc, yc, 1, kfft3, kfft3, 1 } 
           { imin, imax-1, 3, yc-10, yc+10, 3, kfft3, kfft3, 1 }    
      }
   }
}
--- CREATE: config.<part>.in
cfg:CREATE()