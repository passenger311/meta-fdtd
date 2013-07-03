--- Real values in nm

real_length = 100
real_SC = 40
real_wavelength = 1000.0    -- real wavelength

n_bg = 3.42

--- conversion to computation scale

resolution = real_wavelength/10 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

--- PML size parameter
size_pml = 12

--- Padding size parameter
size_pad = 2

--- Computational units
size_tfsf = math.floor(real_length/conv/2+0.5) 
size_SC = math.floor(real_SC/conv/2+0.5)

--- Gaussian envelope of injection field
probe_on = true
pinv_wl = inv_wavelength
pampl = 1e5 --1e7   -- electric field strength in V/m
pwidthl = 0.5    -- width in number of periods
poffsetl = 0    -- offset in number of periods
pattackl = 20 -- attack in number of periods
psustainl = 0   -- sustain in number of periods
pdecayl = 2000    -- decay in number of periods
ppsi = 0

--- BulkSC material parameters
gammap = 1e12
M = 5e-10 --check
egap = 1.424
me = 0.067
mh = 0.46
N0 = 1.5e24
gammanr = 1e8
temp = 300
k_max = 2.5e8
numk = 100
--GaAs transparency density is 1.2e18/cm^3=1.2e24/cm^3

--- diagnostics?
diagpspec=true
samp_fft = 2^13 -- sampling of fft

--- Courant factor
dt = 0.999  -- time step length compared to grid step length (--> Courant stability factor)

ncycles_probe = 2^17-1
ncyc_probe_start = 0 --math.floor((attackl+sustainl+decayl)*resolution/dt+.5)
ncycles = ncyc_probe_start+ncycles_probe -- number of cycles

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("tfsf-domain x dir(grid):                  ", size_tfsf)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
