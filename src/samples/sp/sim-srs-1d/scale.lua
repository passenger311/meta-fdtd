pi = 3.141592653589793116

--- Real values in um

real_wavelength = 1.550    -- real wavelength

real_tfsf_inj = -0.02 -- tfsf_inj point
real_fft = 10 -- fourier point

n_bg = 3.47
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium

mat = 'gold' -- gold, silver 

--- conversion to computation scale

resolution = 60


conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale

real_dx = conv

inv_wavelength = conv/real_wavelength -- inverse wavelength

tfsf_inj = math.floor(real_tfsf_inj/conv+.5)
fft = math.floor(real_fft/conv+.5)

real_c_si = 2.99792458e8
real_eps0_si = 8.8541878e-12 -- []
real_c = 299.792458   --  [um/ps]
 
real_dt = real_dx / real_c; 


--- for chi3 nonlinear material

-- convert frequency [Thz] to inverse wavelength 1/lambda [1/um]


real_invlambdaR  = 15.6 / real_c        -- inverse Raman wavelength [1/um] 
real_gammaR      = math.pi * 0.105      -- Raman damping constant [1/ps]
real_chi3R   = 11.2e-18                 -- Raman susceptibility in [m^2/V^2]

-- 11.2e−14 i cm^2/V^2

conv_e = math.sqrt(real_eps0_si) / ( 1e-6 )^1.5 / real_c_si 


-- Gaussian envelope of injection field
widthl = 150    -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 500 -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 500    -- decay in number of periods
nrefr = n_bg -- reference injection refractive index

ampl_s = 1/100;    -- in units of 1/sqrt(chi3R)
ampl_p = 50/100;

-- PML size parameter
size_pml = 12

-- Padding size parameter
size_pad = 3

-- Courant factor
dt = .95 -- 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 4*32768-1 -- number of cycles

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
