cfg = CONFIG{scenes=true}
--taskid = select(1,...)
cyl_rad = 1
dofile("scale.lua")
dofile("neff.lua")
n_bg = nrefr

hdist_ntff = 5
hdist_tfsf = 5
ncyc_probe_start = 0
ncycles = ncycles_probe
pattackl = pattackl - (attackl+sustainl+decayl)

imin = -1
imax = 1
jmin = -1
jmax = 1
kmin = -hdist_ntff-size_pad
kmax = hdist_ntff+size_pad

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)


--- GRID Definition

cfg:GRID{

   dim = 3,                         -- number of dimensions
   partition = { 0, 1 },
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step length compared to grid step length (--> Courant stability factor)
   irange = { imin, imax },         -- range of computational window in x direction
   jrange = { jmin, jmax },         -- -"- in y direction
   krange = { kmin-size_pml, kmax+size_pml },    -- -"- in z direction
   dx = { conv*1e-9, 1, 1, 1 }

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
   },
   OUT{
      file = { "GPL", "Ex" },
      type = { "Ex", "S", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         BOX{
            { 0, 0, 1, 0, 0, 1, 0, 0, 1 }
         }
      },
--      on = false
   },
   OUT{
      file = { "GPL", "Eap" },
      type = { "Eap", "S", ".F." },
      time = { 0, ncycles, 2000 },
      REG{
         BOX{
            { 0, 0, 1, 0, 0, 1, 0, 0, 1 }
         }
      },
      on = false
   },
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
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      planewave = { phi=0, theta=0, psi=ppsi, nrefr=nrefr },
      config = {0,0,0,0,1,1}
   },
   REG{
      BOX{
         {-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf,hdist_tfsf,1}
      }
   },
   on = probe_on or probe_on_injx
}
cfg:SRC{
   TFSFBOX{
      invlambda = pinv_wl,
      amplitude = pampl*(math.sqrt(8.854187817e-12)*(conv*1e-9)^(3/2)/(frequ_factor*1000)),
      pulse = {
         shape="Gaussian",
         width=math.floor(pwidthl/pinv_wl/dt+.5),
         offset=math.floor(poffsetl/pinv_wl/dt+.5),
         attack=math.floor(pattackl/pinv_wl/dt+.5),
         sustain=math.floor(psustainl/pinv_wl/dt+.5),
         decay=math.floor(pdecayl/pinv_wl/dt+.5)
      },
      planewave = { phi=0, theta=0, psi=ppsi-90, nrefr=nrefr },
      config = {0,0,0,0,1,1}
   },
   REG{
      BOX{
         {-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf-1,hdist_tfsf+1,1,-hdist_tfsf,hdist_tfsf,1}
      }
   },
   on = probe_on_injy
}

cfg:DIAG{
   PSPEC{
      file = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -hdist_tfsf+2, -hdist_tfsf+2, 1 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      mode = "F"
   },
   REG{
      BOX{
         { 0, 0, 1, 0, 0, 1, -hdist_tfsf+2, -hdist_tfsf+2, 1 }
      }
   },
   on = diagmode
}

--- CREATE: config.<part>.in
cfg:CREATE()
