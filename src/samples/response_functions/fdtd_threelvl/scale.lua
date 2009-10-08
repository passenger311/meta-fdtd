pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 500.0    -- real wavelength

real_tfsf_inj = -2 -- tfsf_inj point
real_fft = 102 -- fourier point

n_bg = 1.0
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium
mat = 'gold' -- gold, silver 

--- conversion to computation scale

resolution = 500/1 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

tfsf_inj = math.floor(real_tfsf_inj/conv+.5)
fft = math.floor(real_fft/conv+.5)


-- Gaussian envelope of injection field
widthl = 1    -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 15  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 100    -- decay in number of periods
nrefr = n_bg -- reference injection refractive index

-- PML size parameter
size_pml = 12

-- Padding size parameter
size_pad = 3

-- Courant factor
dt = .9 -- 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 8*32768-1 -- number of cycles

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("tfsf-injection (grid):                    ", tfsf_inj)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)

foutput = io.open("data.save","w+")
foutput:write(real_fft-2,"\n")
foutput:write(conv,"\n")
foutput:write(inv_wavelength,"\n")
foutput:close()
