
if taskid then
   print("taskid = ",taskid)
-- do something like this:
   real_hdnp = 20. + (taskid-1) * 20.-- half distance between 2 nanoparticles

-- for array job set real_length_mmi=scanval.
else
   real_hdnp = 40.-- half distance between 2 nanoparticles

end


pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 500.0    -- real wavelength

real_hhwaveguide = 300.0   -- waveguide half height
real_rwaveguide = 120.0    -- waveguide radius
real_hhcavity = 500.0      -- cavity half height
real_rcavity = 240         -- cavity radius
real_rnp = 20              -- nanoparticle radius

real_kmax = 4*real_hhwaveguide+2*real_hhcavity

real_kinj = 50

real_kfft2 = real_kmax / 2
real_kfft3 = real_kmax - 2*real_hhwaveguide + 100
real_kfft4 = real_kmax - 100

n_waveguide = 3.47   -- waveguide refractive index
n_max = n_waveguide  -- maximum refraxctive index to determine optically thickest medium

-- Drude material in THz
--eps_infD = 9.0685
--real_omegaD = 2155.6     -- Drude plasma frequency [2 pi c]
--real_gammaD = 18.36*2*pi  -- Drude damping constant [1/dt]

-- Drude-Lorentz material in THz
eps_infDL = 5.9673
real_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
real_gammaDL = 15.92*2*pi   -- Drude damping constant [1/dt]
real_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
real_gammaL = 104.86*2*pi   -- Lorentzian damping constant [1/dt]
deltaepsl = 1.09


--- conversion to computation scale

resolution = 40 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hhwaveguide = math.floor(real_hhwaveguide/conv+.5)
rwaveguide = math.floor(real_rwaveguide/conv+.5)
hhcavity = math.floor(real_hhcavity/conv+.5) -- number of grid points representing the cavity half height
rcavity = math.floor(real_rcavity/conv+.5)
rnp = math.floor(real_rnp/conv+.5)
hdnp = math.floor(real_hdnp/conv+.5)

-- Gaussian envelope of injection field
kinj = math.floor(real_kinj/conv+.5)
widthl = 2     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4    -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 4     -- decay in number of periods
nrefr = 3.1378985767173 -- reference injection refractive index

-- Fourier transform parameters
kfft0 = kinj - 2
kfft1 = kinj + 2
kfft2 = math.floor(real_kfft3/conv+.5)
kfft3 = math.floor(real_kfft3/conv+.5)
kfft4 = math.floor(real_kfft4/conv+.5)

-- PML size parameter
size_pml = 11

-- Courant factor
dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 8000 -- number of cycles


field_out = "field" -- geo output file to compute injection profile with Matlab
field_inj = "tfsfInjProfile.set" -- file created by Matlab

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("Half-Height of waveguide in cells (grid): ", hhwaveguide)
print("Radius of waveguide in cells (grid):      ", rwaveguide)
print("Half-Height of cavity in cells (grid):    ", hhcavity)
print("Radius of cavity in cells (grid):         ", rcavity)
print("The following two parameters are only applicable in the case of two nanoparticles")
print("Radius of nanoparticles (grid):           ", rnp)
print("Half-Distance of nanoparticles (grid):    ", hdnp)
