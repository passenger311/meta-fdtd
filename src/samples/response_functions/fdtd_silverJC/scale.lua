pi = 3.141592653589793116

--- Real values in nm

real_wavelength = 500.0    -- real wavelength

real_hdist = 70 -- half distance of capacitor plates
real_rnp = 10

n_bg = 1.0
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium
mat = 'silver' -- gold, silver 

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
end

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
foutput:write(eps_infDL,"\n")
foutput:write(real_omegaDL/frequ_factor*conv,"\n")
foutput:write(real_gammaDL/frequ_factor*conv,"\n")
foutput:write(deltaepsl,"\n")
foutput:write(real_omegaL/frequ_factor*conv,"\n")
foutput:write(real_gammaL/frequ_factor*conv,"\n")
if (mat == 'silver') then
foutput:write(deltaepsl2,"\n")
foutput:write(real_omegaL2/frequ_factor*conv,"\n")
foutput:write(real_gammaL2/frequ_factor*conv,"\n")
end
foutput:close()