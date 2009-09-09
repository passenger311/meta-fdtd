pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 600.0    -- real wavelength

real_hdist = 40.
real_hwidth = 20.

n_bg = 1.0
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium
mat = 'dielectric' -- gold, silver, carbon, dielectric with n_sphere
n_sphere = 6

step_fft = 2

--- conversion to computation scale

resolution = 600/4 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hdist = math.floor(real_hdist/conv+.5)
hwidth = math.floor(real_hwidth/conv+.5)
hdist_tfsf_ij = hwidth
hdist_ntff_ij = hdist_tfsf_ij
hdist_tfsf_k = hdist
hdist_ntff_k = hdist_tfsf_k + 2

-- Gaussian envelope of injection field
widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index
field_inj = "tfsfInjProfile.set" -- file created by Matlab

-- PML size parameter
size_pml = 9

-- Courant factor
dt = .7 -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 8*16384-1 -- number of cycles


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
  real_gammaDL = 12.6126e-12 -- mistake in McMahon et al?
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
print("ntff-domain x/y dir (grid):               ", hdist_ntff_ij)
print("tfsf-domain x/y dir(grid):                ", hdist_tfsf_ij)
print("ntff-domain z dir (grid):                 ", hdist_ntff_k)
print("tfsf-domain z dir (grid):                 ", hdist_tfsf_k)
print("Half length of system (grid):             ", hdist)
print("Half width of system (grid):              ", hwidth)
print("Size of PML (grid):                       ", size_pml)

