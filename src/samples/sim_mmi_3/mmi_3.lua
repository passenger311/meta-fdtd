
cfg = CONFIG{scenes=true}


--- real device parameters (um)

real_wavelength = 1.550

real_width_mmi = 2.43
real_length_mmi = 6.15
real_width_wg = 0.400
real_height_bsio2 = 0.7
real_height_wg = 0.220
real_length_wg1 = 1.4
real_hwidth_sep = 0.65

nref_si = 3.47
nref_sio2 = 1.444

eps_si = nref_si^2
eps_sio2 = nref_sio2^2

--- scaling parameters

resolution = 40

dt = 0.574

real_dx = real_wavelength / nref_si / resolution

print("resolution = ", resolution)
print("real_dx = ", real_dx)


print("wavelength (grid) = ", real_wavelength/real_dx)

invwavelength = real_dx / real_wavelength

print("invwavelength (grid) = ", invwavelength)


--- device in natural HL, dx=dy=dz=1.

ncyc = 8000

pulsehwhm = 500
pulsehsteps = 1000
kpulse = 30

-- grid dimensions: i,j,k

width = 281
height = 121
length = 850

maxi = math.floor( (width - 1) / 2  )
maxj = math.floor( height - 1 )
maxk = math.floor( length - 1 )


height_wg = math.floor(real_height_wg/real_dx+0.5)
height_bsio2 = math.floor(real_height_bsio2/real_dx+0.5)
hwidth_wg = math.floor(real_width_wg/2/real_dx+0.5)
hwidth_mmi = math.floor(real_width_mmi/2/real_dx+0.5)
length_wg1 = math.floor(real_length_wg1/real_dx+0.5)
length_mmi = math.floor(real_length_mmi/real_dx+0.5)
hwidth_sep = math.floor(real_hwidth_sep/real_dx+0.5)
length_wg2 = maxk - length_wg1 - length_mmi


print("irange (grid) = ", -maxi, maxi)
print("jrange (grid) = ", 0, maxj)
print("krange (grid) = ", 0, maxk)
print("height_wg (grid) = ", height_wg)
print("height_bsio2 (grid) = ", height_bsio2)
print("hwidth_wg (grid) = ", hwidth_wg)
print("length_wg1 (grid) = ", length_wg1)
print("length_mmi (grid) = ", length_mmi)
print("hwidth_mmi (grid) = ", hwidth_mmi)
print("hwidth_sep (grid) = ", hwidth_sep)
print("length_wg2 (grid) = ", length_wg2)

print("irange (real) = ", -maxi*real_dx, maxi*real_dx)
print("jrange (real) = ", 0, maxj*real_dx)
print("krange (real) = ", 0, maxk*real_dx)


--- create dielectric structure

scene1 = Scene{value=1.}
box_bsio2 = Box{ 
   from={-maxi-1,-1,-1},
   to={maxi+1,height_bsio2,maxk+1}
}

box_wg1 = Box{
   from={-hwidth_wg,height_bsio2,-1},
   to={hwidth_wg,height_bsio2+height_wg,length_wg1 }
}

box_mmi = Box{
   from={-hwidth_mmi,height_bsio2,length_wg1},
   to={hwidth_mmi,height_bsio2+height_wg,length_wg1+length_mmi }
}

i0 = -hwidth_sep-hwidth_wg
i1 = -hwidth_sep+hwidth_wg

box_wg2 = Box{
   from={i0,height_bsio2,length_wg1+length_mmi},
   to={i1,height_bsio2+height_wg,maxk+1 }
}

i0 = hwidth_sep-hwidth_wg
i1 = hwidth_sep+hwidth_wg

box_wg3 = Box{
   from={i0,height_bsio2,length_wg1+length_mmi},
   to={i1,height_bsio2+height_wg,maxk+1 }
}

scene1:add{ box_bsio2, depth=1, value=eps_sio2 }
scene1:add{ box_wg1, depth=1, value=eps_si }
scene1:add{ box_mmi, depth=1, value=eps_si }
scene1:add{ box_wg2, depth=1, value=eps_si }
scene1:add{ box_wg3, depth=1, value=eps_si }

grid1 = Grid{from={-maxi,height_bsio2-10,0},to={maxi,height_bsio2+height_wg+1,maxk}}
pad = 50
grid2 = Grid{from={-hwidth_wg-pad,height_bsio2-pad,kpulse }, to={hwidth_wg+pad,height_bsio2+height_wg+pad,kpulse}}
grid3 =  Grid{yee=false,from={-maxi,height_bsio2-10,0},to={maxi,height_bsio2+height_wg+1,maxk},offset={-maxi,0,0},cells={50,100,100}}

cfg:CREATE_GEO{"scene1", scene=scene1, grid=grid1, method="default",comps=3, silent=false, on=true }

cfg:CREATE_GEO{"inj1", scene=scene1, grid=grid2, method="default",comps=3, silent=false, on=true }

cfg:CREATE_PREVIEW{"scene1", scene=scene1, grid=grid3, method="default", silent=false, on=true }

--- GRID Definition

cfg:GRID{

   dim = 3,
   partition = { 0, 1 },
   ncyc = ncyc,
   dt = dt,
--   irange = { -maxi, maxi },
   irange = { 0, maxi },
   jrange = { 0, maxj },
   krange = { 0, maxk }

}
--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{
	    { -maxi-1, maxi+1, 1, 0, height_bsio2-1, 1, 0, maxk, 1, ":", eps_sio2, eps_sio2, eps_sio2 },
	    { -maxi-1, maxi+1, 1, height_bsio2, maxj, 1, 0, maxk, 1, ":", 1.,1.,1. }
	 },
	 LOAD_GEO{ "scene1" }
      },
      on = true
   },

   OUT{
      file = { "VTK", "slice0_xy_eps" },
      type = { "Eps", "N" },
      time = { 0, ncyc, 1000 },
      REG{
	 BOX{
	    { -maxi,maxi, 1, 
	       0, maxj, 1, 
	       31, 31, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "slice0_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       31, 31, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "slice1_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	      math.floor(length_wg1/2), math.floor(length_wg1/2) , 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "slice2_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-10,hwidth_wg+10, 1, 
	       height_bsio2-10, height_bsio2+height_wg+10, 1, 
	       length_wg1, length_wg1, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "slice3_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_mmi-20,hwidth_mmi+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	      length_wg1+math.floor(length_mmi/2),length_wg1+math.floor(length_mmi/2) , 1  }
	 }
      }
   },

    OUT{
       file = { "VTK", "slice4_xy_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
       REG{
	  BOX{
	     { -hwidth_sep-hwidth_wg-20,hwidth_sep+hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       length_wg1+length_mmi,length_wg1+length_mmi, 1  }
	  }
       }
    },
    
    OUT{
       file = { "VTK", "slice5_xy_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
       REG{
	  BOX{
	     { -hwidth_sep-hwidth_wg-20,hwidth_sep+hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       length_wg1+length_mmi+math.floor(length_wg2/2),length_wg1+length_mmi+math.floor(length_wg2/2) , 1  }
	  }
       }
    },

      
    OUT{
       file = { "VTK", "slice1_xz_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
       REG{
	  BOX{
	     { -maxi+11, maxi-11, 1, 
	       height_bsio2+math.floor(height_wg/2),height_bsio2+math.floor(height_wg/2), 1,
	       0+11, maxk-11, 1	       
	    }
	  }
       }
    },

}

--- BOUND Definition

cfg:BOUND{
   config = { 2, 1, 1, 1, 1, 1 },
   PML{
      cells = 11,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFINJ{
      invlambda = invwavelength,
      amplitude = 1.,
      pulse = { 
	 shape="Gaussian", 
	 width=pulsehwhm,
	 offset=0, attack=pulsehsteps, sustain=0, decay=pulsehsteps
      },
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=nref_si }
   },
   REG{
      LOAD{ "inj1.set" }
   },
   on = true
}


--- PSPEC

yc = height_bsio2+math.floor(height_wg/2)

cfg:DIAG{
   PSPEC{
      file = "input",
      time = { 1, 8192, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { 1, 21, 3, yc-10, yc+10, 3, 90, 90, 1 }  
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "output",
      time = { 1, 8192, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, 700, 700, 1 }  
      }
   }
}


--- CREATE: config.<part>.in
cfg:CREATE()