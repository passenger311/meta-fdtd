
cfg = CONFIG{scenes=true}
dofile("scale.lua")

hdist_ntff = 10
hdist_tfsf = 10

imin = -hdist_ntff
imax = hdist_ntff
jmin = -hdist_ntff
jmax = hdist_ntff
kmin = -hdist_ntff-size_pad
kmax = hdist_ntff+size_pad

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(n_max*math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 3,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin, imax },   -- range of computational window in x direction
   jrange = { jmin, jmax },   -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml }    -- -"- in z direction

}

--- FDTD Definition
eps_bg = n_bg^2
cfg:FDTD{

   EPSILON{
      REG{
         BOX{
            { imin-size_pml-1, imax+size_pml+1, 1, jmin-size_pml-1, jmax+size_pml+1, 1, kmin-size_pml-1, kmax+size_pml+1, 1, ":", eps_bg, eps_bg, eps_bg }
         }
      },
      on = true
   }

}

--- BOUND Definition

cfg:BOUND{

   config = { 3, 3, 3, 3, 1, 1 },

   PML{
      cells = size_pml,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFBOX{
      invlambda = inv_wavelength,
      amplitude = 1.0,
      pulse = { 
         shape="Gaussian",
         width=math.floor(widthl*resolution*n_max+.5),
         offset=math.floor(offsetl*resolution*n_max+.5),
         attack=math.floor(attackl*resolution*n_max+.5),
         sustain=math.floor(sustainl*resolution*n_max+.5),
         decay=math.floor(decayl*resolution*n_max+.5)
      },
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=nrefr },
      config = { 0, 0, 0, 0, 1, 1}
   },
   REG{
      BOX{
         {-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf,hdist_tfsf,1}
      }
   },
   on = true
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref-zabs",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf+2, hdist_tfsf-2, 1, -hdist_tfsf+2, hdist_tfsf-2, 1, -hdist_tfsf+2, -hdist_tfsf+2, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "EHT",
      time = { 0, ncycles, (ncycles+1)/4096 },
      mode = "EHT"
   },
   REG{
      BOX{
         { -1, 1 , 1, -1, 1, 1, 0, 0, 1 }
      }
   }
}

--- CREATE: config.<part>.in
cfg:CREATE()
