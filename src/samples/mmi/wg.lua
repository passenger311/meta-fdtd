
cfg = CONFIG{scenes=true}


--- real device parameters (um)

real_wavelength = 1.550

real_width_wg = 0.400
real_height_bsio2 = 0.7
real_height_wg = 0.220
real_kinj = 0.2
real_kfft1 = real_kinj + 1.6 
real_kfft2 = - 0.4

nref_si = 3.47
nref_sio2 = 1.444

eps_si = nref_si^2
eps_sio2 = nref_sio2^2

--- scaling parameters

resolution = 35

dt = 0.574

real_dx = real_wavelength / nref_si / resolution

print("resolution = ", resolution)
print("real_dx = ", real_dx)


print("wavelength (grid) = ", real_wavelength/real_dx)

invwavelength = real_dx / real_wavelength

print("invwavelength (grid) = ", invwavelength)


--- device in natural HL, dx=dy=dz=1.

ncyc = 16384

pulsehwhm = 150
pulsehsteps = 500
kpulse = 30

-- grid dimensions: i,j,k

width = 101
height = 121
length = 850

imax = math.floor( (width - 1) / 2  )
jmax = math.floor( height - 1 )
kmax = math.floor( length - 1 )


height_wg = math.floor(real_height_wg/real_dx+0.5)
height_bsio2 = math.floor(real_height_bsio2/real_dx+0.5)
hwidth_wg = math.floor(real_width_wg/2/real_dx+0.5)
length_wg1 = length
hwidth_sep = 0
kinj =  math.floor(real_kinj/real_dx+0.5)
kfft1 =  math.floor(real_kfft1/real_dx+0.5)
kfft2 =  kmax + math.floor(real_kfft2/real_dx+0.5)
yc = height_bsio2+math.floor(height_wg/2+0.5)


print("irange (grid) = ", -imax, imax)
print("jrange (grid) = ", 0, jmax)
print("krange (grid) = ", 0, kmax)
print("height_wg (grid) = ", height_wg)
print("height_bsio2 (grid) = ", height_bsio2)
print("hwidth_wg (grid) = ", hwidth_wg)
print("length_wg1 (grid) = ", length_wg1)
print("kinj (grid) = ", kinj)
print("kfft1 (grid) = ", kfft1)
print("kfft2 (grid) = ", kfft2)
print("yc (grid) = ", yc)

print("irange (real) = ", -imax*real_dx, imax*real_dx)
print("jrange (real) = ", 0, jmax*real_dx)
print("krange (real) = ", 0, kmax*real_dx)


--- create dielectric structure

wg = Scene{value=1.}
box_bsio2 = Box{ 
   from={-imax-1,-1,-1},
   to={imax+1,height_bsio2,kmax+1}
}

box_wg1 = Box{
   from={-hwidth_wg,height_bsio2,-1},
   to={hwidth_wg,height_bsio2+height_wg,length_wg1+1 }
}

wg:add{ box_bsio2, depth=1, value=eps_sio2 }
wg:add{ box_wg1, depth=1, value=eps_si }

grid_eps = Grid{from={-imax,height_bsio2-10,0},to={imax,height_bsio2+height_wg+1,kmax}}
pad = 50
grid_inj = Grid{from={-hwidth_wg-pad,height_bsio2-pad,kpulse }, to={hwidth_wg+pad,height_bsio2+height_wg+pad,kpulse}}
grid_prev =  Grid{yee=false,from={-imax,height_bsio2-10,0},to={imax,height_bsio2+height_wg+1,kmax},offset={-imax,0,0},cells={50,100,100}}

cfg:CREATE_GEO{"wg", scene=wg, grid=grid_eps, method="default",comps=3, silent=false, on=true }
cfg:CREATE_GEO{"inj", scene=wg, grid=grid_inj, method="default",comps=3, silent=false, on=true }
cfg:CREATE_PREVIEW{"wg", scene=wg, grid=grid_prev, method="default", silent=false, on=true }

--- GRID Definition

cfg:GRID{

   dim = 3,
   partition = { 0, 1 },
   ncyc = ncyc,
   dt = dt,
--   irange = { -imax, imax },
   irange = { 0, imax },
   jrange = { 0, jmax },
   krange = { 0, kmax }

}
--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{
	    { -imax-1, imax+1, 1, 0, height_bsio2-1, 1, 0, kmax, 1, ":", eps_sio2, eps_sio2, eps_sio2 },
	    { -imax-1, imax+1, 1, height_bsio2, jmax, 1, 0, kmax, 1, ":", 1.,1.,1. }
	 },
	 LOAD_GEO{ "wg" }
      },
      on = true
   },

   OUT{
      file = { "GPL", "wg_point_e_input" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 0, yc, kinj+1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "wg_point_e_middle" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 0, yc, math.floor(kmax/2) }  
	 }
      }
   },

   OUT{
      file = { "GPL", "wg_point_e_output" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { hwidth_sep, yc, kmax-11 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "wg_point_en_fft1" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { 1, 21, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "wg_point_en_fft2" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }  
	 }
      }
   },

   OUT{
      file = { "VTK", "wg_slice0_xy_eps" },
      type = { "Eps", "N" },
      time = { 0, ncyc, 1000 },
      REG{
	 BOX{
	    { -imax,imax, 1, 
	       0, jmax, 1, 
	       kinj+1, kinj+1, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "wg_slice0_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       kinj-1, kinj-1, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "wg_slice1_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       kinj+1, kinj+1, 1  }
	 }
      }
   },

   OUT{
      file = { "VTK", "wg_slice2_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	      kfft1, kfft1, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "wg_slice3_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
      REG{
	 BOX{
	    { -hwidth_wg-20,hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	      math.floor(length_wg1/2),math.floor(length_wg1/2) , 1  }
	 }
      }
   },

    OUT{
       file = { "VTK", "wg_slice4_xy_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
       REG{
	  BOX{
	     { -hwidth_sep-hwidth_wg-20,hwidth_sep+hwidth_wg+20, 1, 
	       height_bsio2-20, height_bsio2+height_wg+20, 1, 
	       kfft2,kfft2, 1  }
	  }
       }
    },
      
    OUT{
       file = { "VTK", "wg_slice1_xz_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 1000 },
       REG{
	  BOX{
	     { -imax+11, imax-11, 1, 
	       yc, yc, 1,
	       0+11, kmax-11, 1	       
	    }
	  }
       }
    }

}

--- BOUND Definition

cfg:BOUND{
   config = { 0, 1, 1, 1, 1, 1 },
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
      LOAD{ "tfsfinj.set" }
   },
   on = true
}


--- PSPEC


cfg:DIAG{
   PSPEC{
      file = "wg_in",
      time = { 1, 16384, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { 1, 21, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 }  
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "wg_out",
      time = { 1, 16384, 1 },
      reffile = "wg_out",
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }  
      }
   }
}


--- CREATE: config.<part>.in
cfg:CREATE()