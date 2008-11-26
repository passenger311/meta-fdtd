
cfg = CONFIG{scenes=true}

taskid = select(1,...)

dofile("scale.lua")

-- grid dimensions: i,j,k

width = math.floor((2*real_width_wg+real_width_mmi)/real_dx)

imax = math.floor( (width - 1) / 2  )
jmax = math.floor( height - 1 )
kmax = math.floor( length - 1 )

print("irange (grid) = ", -imax, imax)
print("jrange (grid) = ", 0, jmax)
print("krange (grid) = ", 0, kmax)

print("irange (real) = ", -imax*real_dx, imax*real_dx)
print("jrange (real) = ", 0, jmax*real_dx)
print("krange (real) = ", 0, kmax*real_dx)



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


--- create dielectric structure

mmi = Scene{ value=1.}


wg = Scene{value=1.}
box_bsio2 = Box{ 
   from={-imax0-100,jmin0-100,kmin0},
   to={imax0+100,height_bsio2,kmax0+1}
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
   to={i1,height_bsio2+height_wg,kmax+1 }
}

i0 = hwidth_sep-hwidth_wg
i1 = hwidth_sep+hwidth_wg

box_wg3 = Box{
   from={i0,height_bsio2,length_wg1+length_mmi},
   to={i1,height_bsio2+height_wg,kmax+1 }
}

mmi:add{ box_bsio2, depth=1, value=eps_sio2 }
mmi:add{ box_wg1, depth=1, value=eps_si }
mmi:add{ box_mmi, depth=1, value=eps_si }
mmi:add{ box_wg2, depth=1, value=eps_si }
mmi:add{ box_wg3, depth=1, value=eps_si }

grid_eps = Grid{from={-imax,height_bsio2-10,0},to={imax,height_bsio2+height_wg+1,kmax}}
pad = 50
grid_inj = Grid{from={-hwidth_wg-pad,height_bsio2-pad,kinj }, to={hwidth_wg+pad,height_bsio2+height_wg+pad,kinj}}
grid_prev =  Grid{yee=false,from={-imax,height_bsio2-10,0},to={imax,height_bsio2+height_wg+1,kmax},offset={-imax,0,0},cells={50,100,100}}

cfg:CREATE_GEO{"mmi", scene=mmi, grid=grid_eps, method="default",comps=3, silent=false, on=true }
cfg:CREATE_GEO{"inj", scene=mmi, grid=grid_inj, method="default",comps=3, silent=false, on=true }
cfg:CREATE_PREVIEW{"mmi", scene=mmi, grid=grid_prev, method="default", silent=false, on=true }

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
	    { imin0, imax0+1, 1, jmin0, height_bsio2-1, 1, kmin0, kmax0+1, 1, ":", eps_sio2, eps_sio2, eps_sio2 },
	    { imin0, imax0+1, 1, height_bsio2, jmax0+1, 1, kmin0, kmax0+1, 1, ":", 1.,1.,1. }
	 },
	 LOAD_GEO{ "mmi" }
      },
      on = true
   },


   OUT{
      file = { "SET", "mmi_eps" },
      type = { "Eps", "N" },
      time = { 0, 0, 1 },
      REG{
	 BOX{
	    {  imin,imax, 1, 
	       jmin,jmax, 1, 
	       kinj, kinj, 1  }
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_injb" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 0, yc, kinj-1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_injf" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 0, yc, kinj+1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_mid" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { 0, yc, math.floor(kmax/2) }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_e_end" },
      type = { "E", "N", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 POINT{ 
	    { hwidth_sep, yc, kmax }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft0" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { 1, 21, 3, yc-10, yc+10, 3, kinj+1, kinj+1, 1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft1" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { 1, 21, 3, yc-10, yc+10, 3, kfft1, kfft1, 1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft2" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }  
	 }
      }
   },

   OUT{
      file = { "GPL", "mmi_point_en_fft3" },
      type = { "En", "S", ".F." },
      time = { 0, ncyc, 10 },
      REG{
	 BOX{ 
	    { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft3, kfft3, 1 }  
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
      time = { 0, ncyc, 500 },
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
      time = { 0, ncyc, 500  },
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
      time = { 1000, ncyc, 500 },
      REG{
	 BOX{
	    {  imin, imax, 1, 
	       jmin, jmax, 1, 
	      kfft1, kfft1, 1  }
	 }
      }
   },

    OUT{
      file = { "VTK", "mmi_slice3_xy_e" },
      type = { "E", "N" },
      time = { 1000, ncyc, 500 },
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
      time = { 1000, ncyc, 500 },
       REG{
	  BOX{
	     { imin, imax, 1, 
	       jmin, jmax, 1, 
	       kfft2,kfft2, 1  }
	  }
       }
    },

    OUT{
       file = { "VTK", "mmi_slice5_xy_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 500 },
       REG{
	  BOX{
	     { imin, imax, 1, 
	       jmin, jmax, 1, 
	       kfft3,kfft3, 1  }
	  }
       }
    },
      
    OUT{
       file = { "VTK", "mmi_slice1_xz_e" },
       type = { "E", "N" },
      time = { 1000, ncyc, 500 },
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
      file = "mmi_fft0",
      time = { 1, ncyc, 1 },
      mode = "Eap",
      phasewrap = { 1, 0 },
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { 1, 21, 3, yc-10, yc+10, 3, kinj+1, kinj+1, 1 }  
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "mmi_fft1",
      time = { 1, ncyc, 1 },
      mode = "Eap",
      phasewrap = { 1, 0 },
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
      file = "mmi_fft2",
      time = { 1, ncyc, 1 },
--      reffile = "wg3",
      phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft2, kfft2, 1 }  
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "mmi_fft3",
      time = { 1, ncyc, 1 },
--      reffile = "wg2",
       phasewrap = { 1, 0 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{ 
	 { hwidth_sep+1, hwidth_sep+21, 3, yc-10, yc+10, 3, kfft3, kfft3, 1 }  
      }
   }
}


--- CREATE: config.<part>.in
cfg:CREATE()