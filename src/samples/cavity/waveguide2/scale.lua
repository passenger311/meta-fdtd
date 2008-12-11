if (taskid and (taskid ~= "undefined")) then
   print("taskid = ",taskid)
-- do something like this:
--   real_hdnp = 20. + (taskid-1) * 20.-- half distance between 2 nanoparticles
  dxnp=(taskid-1)*10
-- for array job set real_length_mmi=scanval.
else
-- real_hdnp = 40.-- half distance between 2 nanoparticles
  dxnp=30
end

pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 500.0    -- real wavelength

real_hhwaveguide = 800.0   -- waveguide half height
real_rwaveguide = 120.0    -- waveguide radius
real_rnp = {20,20,20,20,20,20,20,20}        -- nanoparticle radius
real_xnp = {20,0,-20,0,14,14,-14,-14}        -- x position of center of nanoparticle relative to wg center
real_ynp = {0,20,0,-20,14,-14,-14,14}        -- y position of center of nanoparticle
real_znp = {20,20,20,20,-20,-20,-20,-20}        -- z position of center of nanoparticle

real_kmax = 2*real_hhwaveguide

real_kinj = 100


real_kfft0 = real_kinj - 10
real_kfft1 = real_kinj + 10
real_kfft2 = real_kmax * 1/4
real_kfft3 = real_kmax / 2
real_kfft4 = real_kmax *3/4
real_kfft5 = real_kmax - real_kinj

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

resolution = 50 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hhwaveguide = math.floor(real_hhwaveguide/conv+.5)
rwaveguide = math.floor(real_rwaveguide/conv+.5)
rnp={}; inp={}; jnp={}; knp={}
for i,v in ipairs(real_rnp) do
  rnp[i] = math.floor(real_rnp[i]/conv+.5)
  inp[i] = math.floor(real_xnp[i]/conv+.5)
  jnp[i] = math.floor(real_ynp[i]/conv+.5)
  knp[i] = math.floor(real_znp[i]/conv+.5)
end

-- Gaussian envelope of injection field
kinj = math.floor(real_kinj/conv+.5)
widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 10     -- decay in number of periods
nrefr = 1 --3.1493530190359 -- reference injection refractive index

-- Fourier transform parameters
kfft0 = math.floor(real_kfft0/conv+.5)
kfft1 = math.floor(real_kfft1/conv+.5)
kfft2 = math.floor(real_kfft2/conv+.5)
kfft3 = math.floor(real_kfft3/conv+.5)
kfft4 = math.floor(real_kfft4/conv+.5)
kfft5 = math.floor(real_kfft5/conv+.5)

-- PML size parameter
size_pml = 11

-- Padding size parameter
size_pad = 20

-- Courant factor
dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 10000 -- number of cycles


field_out = "field" -- geo output file to compute injection profile with Matlab
field_inj = "tfsfInjProfile.set" -- file created by Matlab

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("Half-Height of waveguide (grid):          ", hhwaveguide)
print("Radius of waveguide(grid):                ", rwaveguide)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
print("The following parameters are only applicable in the case of nanoparticles")
for i,v in ipairs(rnp) do
print("Radius of nanoparticle (grid):            ", v)
print("Position of nanoparticle (grid):          ", inp[i], jnp[i], knp[i])
end
print("\nPosition of Injection plane (grid):       ", kinj)
print("Position of Fourier transforms (grid):    ", kfft0, kfft1, kfft2, kfft3, kfft4, kfft5)
