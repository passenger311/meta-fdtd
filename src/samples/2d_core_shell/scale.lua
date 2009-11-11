if (taskid and (taskid ~= "undefined")) then
   print("taskid = ",taskid)
else
   taskid = 1 
end

pi = 3.141592653589793116


--- Real values in nm

real_wavelength = 650.0    -- real wavelength

real_core = {40} -- total radius of shell particle
real_shell = {6} -- size of shell
real_xnp = {0} 
real_ynp = {0} 
real_znp = {0}

real_hdist_tfsf_i = 48
real_hdist_tfsf_j = 48
real_hdist_tfsf_k = 2

n_diel_core = math.sqrt(2.04) --silica
n_diel_shell = math.sqrt(2.04)
n_bg = 1.00
n_max = n_bg  -- maximum refractive index to determine optically thickest medium
metal_core = false
metal_shell = false
if (metal_core and metal_shell) then
  mat_core = 'gold' -- gold,silver, carbon
  mat_shell = 'gold' -- gold,silver, carbon
elseif (metal_core and not metal_shell) then
  mat_core = 'gold' -- gold,silver, carbon
  mat_shell = 'diel'
elseif (metal_shell and not metal_core) then
  mat_core = 'diel'
  mat_shell = 'gold' -- gold, silver
else
  mat = 'dielectric'
end

step_dft = 2
sampl_dft = 512
step_fft = 2
sampl_fft = 512


--- conversion to computation scale

resolution = real_wavelength/.5 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hdist_tfsf_i = math.floor(real_hdist_tfsf_i/conv+.5)
hdist_ntff_i = hdist_tfsf_i + 2 --math.floor(real_hdist_ntff_i/conv+.5)
hdist_tfsf_j = math.floor(real_hdist_tfsf_j/conv+.5)
hdist_ntff_j = hdist_tfsf_j + 2 --math.floor(real_hdist_ntff_j/conv+.5)
hdist_tfsf_k = math.floor(real_hdist_tfsf_k/conv+.5)
hdist_ntff_k = hdist_tfsf_k + 2 --math.floor(real_hdist_ntff_k/conv+.5)

rnp={}; shell={}; inp={}; jnp={}; knp={}
for i,v in ipairs(real_core) do
  rnp[i] = math.floor(real_core[i]/conv+.5)
  shell[i] = math.floor(real_shell[i]/conv+.5)
  inp[i] = math.floor(real_xnp[i]/conv+.5)
  jnp[i] = math.floor(real_ynp[i]/conv+.5)
  knp[i] = math.floor(real_znp[i]/conv+.5)
end


--- Gaussian envelope of injection field
widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index
phi=90.0
theta=90.0
psi=90.0


--- PML size parameter

size_pml = 9


--- Padding size parameter

size_pad = 3


--- Courant factor

dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 2*32768-1 -- number of cycles


--- Drude-Lorentz material in THz

if (mat_core == 'gold') then
  eps_infDL = 5.9673
  real_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
  real_gammaDL = 15.92*2*pi   -- Drude damping constant [1/dt]
  real_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
  real_gammaL = 104.86/2*2*pi   -- Lorentzian damping constant [1/dt]
  deltaepsl = 1.09
elseif (mat_core == 'silver') then
  eps_infDL = 1.17152
  real_omegaDL = 13960.4/2/pi
  real_gammaDL = 12.6126e-12 -- mistake in McMahon et al?
  real_omegaL = 8257.18/2/pi
  real_gammaL = 195.614
  deltaepsl = 2.23994
  real_omegaL2 = 3057.07/2/pi 
  real_gammaL2 = 852.675
  deltaepsl2 = 0.222651
elseif (mat_core == 'carbon') then
  eps_infDL = 2.554
  real_omegaDL = 2149.77/2/pi
  real_gammaDL = 9475.68
  real_omegaL = 6820.0/2/pi
  real_gammaL = 4799.37
  deltaepsl = 3.56779
  real_omegaL2 = 2833.44/2/pi 
  real_gammaL2 = 1826.16
  deltaepsl2 = 2.51021
end
if (mat_shell == 'gold') then
  eps_s_infDL = 5.9673
  real_s_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
  real_s_gammaDL = 15.92*2*pi   -- Drude damping constant [1/dt]
  real_s_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
  real_s_gammaL = 104.86/2*2*pi   -- Lorentzian damping constant [1/dt]
  s_deltaepsl = 1.09
elseif (mat_shell == 'silver') then
  eps_s_infDL = 1.17152
  real_s_omegaDL = 13960.4/2/pi
  real_s_gammaDL = 12.6126e-12 -- mistake in McMahon et al?
  real_s_omegaL = 8257.18/2/pi
  real_s_gammaL = 195.614
  s_deltaepsl = 2.23994
  real_s_omegaL2 = 3057.07/2/pi
  real_s_gammaL2 = 852.675
  s_deltaepsl2 = 0.222651
elseif (mat_shell == 'carbon') then
  eps_s_infDL = 2.554
  real_s_omegaDL = 2149.77/2/pi
  real_s_gammaDL = 9475.68
  real_s_omegaL = 6820.0/2/pi
  real_s_gammaL = 4799.37
  s_deltaepsl = 3.56779
  real_s_omegaL2 = 2833.44/2/pi
  real_s_gammaL2 = 1826.16
  s_deltaepsl2 = 2.51021
end

eps_diel_core = n_diel_core^2
eps_diel_shell = n_diel_shell^2

--- print some parameters

print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("ntff-domain x dir (grid):                 ", hdist_ntff_i)
print("tfsf-domain x dir (grid):                 ", hdist_tfsf_i)
print("ntff-domain y dir (grid):                 ", hdist_ntff_j)
print("tfsf-domain y dir (grid):                 ", hdist_tfsf_j)
print("ntff-domain z dir (grid):                 ", hdist_ntff_k)
print("tfsf-domain z dir (grid):                 ", hdist_tfsf_k)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
for i,v in ipairs(rnp) do
  print("Radius of nanoparticle (grid):            ", v)
  print("Radius of core of nanoparticle (grid):    ", v-shell[i])
  print("Position of nanoparticle (grid):          ", inp[i], jnp[i], knp[i])
end


--- Write invlambda.in file for DFT

foutput = io.open("invlambda.in","w+")
foutput2 = io.open("lambda.in","w+")
for i = 500,900,5 do
   foutput:write(conv/i, "\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("invlambda2.in","w+")
foutput2 = io.open("lambda2.in","w+")
for i = 500,900,10 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(2*(hdist_ntff_i+2), " ", 2*(hdist_ntff_j+2), " ", 2*(hdist_ntff_k+2), "\n")
foutput:write(2*(hdist_tfsf_i-2), " ", 2*(hdist_tfsf_j-2), " ", 2*(hdist_tfsf_k-2), "\n")
foutput:write(theta, " ", phi, " ", psi)
foutput:close()
