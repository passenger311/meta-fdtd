pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 500.0    -- real wavelength

real_hdist = 70 -- half distance of capacitor plates
real_rnp = 10

n_bg = 1.0
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium
mat = 'gold' -- gold, silver 

--- conversion to computation scale

resolution = 500/1 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hdist = math.floor(real_hdist/conv+.5)
hdist_tfsf_k = hdist + 12

rnp = math.floor(real_rnp/conv+.5)

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

ncycles = 1*32768-1 -- number of cycles


--- Drude-Lorentz material in THz

eps_C_infDL = 2.554
real_C_omegaDL = 2149.77/2/pi
real_C_gammaDL = 9475.68
real_C_omegaL = 6820.0/2/pi
real_C_gammaL = 4799.37
C_deltaepsl = 3.56779
real_C_omegaL2 = 2833.44/2/pi 
real_C_gammaL2 = 1826.16
C_deltaepsl2 = 2.51021

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("tfsf-domain (grid):                       ", hdist_tfsf_k)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(inv_wavelength,"\n")
foutput:write(eps_C_infDL,"\n")
foutput:write(real_C_omegaDL/frequ_factor*conv,"\n")
foutput:write(real_C_gammaDL/frequ_factor*conv,"\n")
foutput:write(C_deltaepsl,"\n")
foutput:write(real_C_omegaL/frequ_factor*conv,"\n")
foutput:write(real_C_gammaL/frequ_factor*conv,"\n")
foutput:write(C_deltaepsl2,"\n")
foutput:write(real_C_omegaL2/frequ_factor*conv,"\n")
foutput:write(real_C_gammaL2/frequ_factor*conv,"\n")
foutput:close()
