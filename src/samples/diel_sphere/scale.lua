pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 700.0    -- real wavelength

real_rsphere = 1500
real_htfsf = 3500

n_sphere = 1.1
n_bg = 1.0
n_max = n_sphere  -- maximum refraxctive index to determine optically thickest medium

--- conversion to computation scale

resolution = 20 -- resolution of wavelength in optically thickest medium

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

htfsf = math.floor(real_htfsf/conv+.5)
rsphere = math.floor(real_rsphere/conv+.5)

--- Gaussian envelope of injection field

widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
nrefr = n_bg --3.1493530190359 -- reference injection refractive index
field_inj = "tfsfInjProfile.set" -- file created by Matlab

--- PML size parameter

size_pml = 9

--- Padding size parameter

size_pad = 5 -- distance of PML to tfsf box

--- Courant factor

dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 4095 -- number of cycles


--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("tfsf-domain (grid):                       ", htfsf)
print("Radius of dielectric sphere (grid):       ", rsphere)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
