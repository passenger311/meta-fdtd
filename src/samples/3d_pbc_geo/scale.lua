pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 550.0    -- real wavelength

real_hxhole = 50
real_hyhole = 55
real_hperiod = 120
real_hholeheight = 20
real_hsio2height = 70
real_hmetalheight = 60

n_bg = 1.0
n_max = n_bg
eps_sio2 = 2.1
mat = 'gold' -- gold, silver, carbon, dielectric with n_sphere

step_fft = 4
samp_fft = 512
samp_dft = 512

--- conversion to computation scale

resolution = 550/2-- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hxhole = math.floor(real_hxhole/conv+.5)
hyhole = math.floor(real_hyhole/conv+.5)

hperiod = math.floor(real_hperiod/conv+.5)
hholeheight = math.floor(real_hholeheight/conv+.5)
hsio2height = math.floor(real_hsio2height/conv+.5)
hmetalheight = math.floor(real_hmetalheight/conv+.5)

hdist_tfsf_i = 2*hperiod -- math.floor(real_hdist_tfsf/conv+.5)
hdist_ntff_i = hdist_tfsf_i -- math.floor(real_hdist_ntff/conv+.5)
hdist_tfsf_j = 2*hperiod -- math.floor(real_hdist_tfsf/conv+.5)
hdist_ntff_j = hdist_tfsf_j -- math.floor(real_hdist_ntff/conv+.5)
hdist_tfsf_k = hholeheight+hsio2height+hmetalheight + 30
hdist_ntff_k = hdist_tfsf_k + 2

--- Gaussian envelope of injection field
widthl = 2     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 8  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 15     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index

--- PML size parameter
size_pml = 12

--- Padding size parameter
size_pad = 20

--- Courant factor
dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 4*16384-1 -- number of cycles


--- Drude-Lorentz material in THz

if (mat == 'gold') then
  eps_infDL = 5.9673
  real_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
  real_gammaDL = 15.92*2*pi   -- Drude damping constant [1/dt]
  real_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
  real_gammaL = 104.86/2*2*pi   -- Lorentzian damping constant [1/dt]
  deltaepsl = 1.09
elseif (mat == 'silver') then
  eps_infDL = 1.17152
  real_omegaDL = 13960.4/2/pi
  real_gammaDL = 12.6126e-12
  real_omegaL = 8257.18/2/pi
  real_gammaL = 195.614
  deltaepsl = 2.23994
  real_omegaL2 = 3057.07/2/pi 
  real_gammaL2 = 852.675
  deltaepsl2 = 0.222651
elseif (mat == 'carbon') then
  eps_infDL = 2.554
  real_omegaDL = 2149.77/2/pi
  real_gammaDL = 9475.68
  real_omegaL = 6654.7/2/pi
  real_gammaL = 425.39
  deltaepsl = 3.56779
  real_omegaL2 = 7874.37/2/pi 
  real_gammaL2 = 832.863
  deltaepsl2 = 2.51021
elseif (mat == 'dielectric') then
  eps_infDL = n_sphere^2
end

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("ntff-domain x dir (grid):                 ", hdist_ntff_i)
print("tfsf-domain x dir(grid):                  ", hdist_tfsf_i)
print("ntff-domain y dir (grid):                 ", hdist_ntff_j)
print("tfsf-domain y dir(grid):                  ", hdist_tfsf_j)
print("ntff-domain z dir (grid):                 ", hdist_ntff_k)
print("tfsf-domain z dir (grid):                 ", hdist_tfsf_k)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)

--- Write invlambda.in file for DFT

foutput = io.open("invlambda2.in","w+")
foutput2 = io.open("lambda2.in","w+")
for i = 550,650,10 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:close()
