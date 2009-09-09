pi = 3.141592653589793116

-- Real dimensions in nm (if dimensions not in nm change frequ_factor as well)

real_wavelength = 500.0 
resolution = 150 -- resolution of wavelength in optically thickest medium

n_cyl = 8         -- number of cylinders
real_dcyl = 40    -- diameter of cylinders
real_sx = 100     -- distance between centers of cylinders in x-direction
real_length = 500      -- total length of structure in y-direction

real_offset1 = 70 -- offset between PML and center of first cylinder
real_offsetN = 70 -- offset between center of last cylinder and PML

inj_field = "TE"  -- "TM" or "TE" for transverse magnetic or transverse electric incoupled fields
real_yinj = 20    -- distance of injection plane from PML (should be smaller than real_offset1-real_dcyl/2)

real_yfft1 = 25   -- distance of first Fourier Transform plane from near PML (should be smaller than real_offset1-real_dcyl/2 and larger than real_yinj)
real_yfft2 = 25   -- distance of second Fourier Transform plane from far PML (should be smaller than real_offset2-real_dcyl/2)

dim = 2           -- number of dimensions

n_bg = 1.41       -- background refractive index
n_max = n_bg      -- maximum refractive index of optically thickest medium
nrefr = n_bg      -- refractive index at injection plane

size_pml = 11     -- PML size parameter

dt = 0.7 -- 0.574 -- time step length compared to grid step length (--> Courant stability factor), has to be changed according to dimension

ncycles = 4095    -- number of cycles of total run

-- the following two parameters are only applicable in three dimensional simulations
real_sz = 130     -- distance between centers of cylinders in z-direction
real_hhcyl = 60   -- half the height of each cylinder (should be smaller than real_sz/2)

-- Gaussian envelope of injection field
widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4    -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
if (inj_field == "TE") then
  psi = 0
elseif (inj_field == "TM") then
  psi = 90
else
  printf("ERROR: Injection field not defined properly")
end

-- Drude material parameters in THz
--eps_infD = 9.0685
--real_omegaD = 2155.6      -- Drude plasma frequency [2 pi c]
--real_gammaD = 18.36*2*pi  -- Drude damping constant [1/dt]

-- Drude-Lorentz material parameters in THz
eps_infDL = 5.9673
real_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
real_gammaDL = 15.92*2*pi  -- Drude damping constant [1/dt]
real_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
real_gammaL = 104.86*2*pi  -- Lorentzian damping constant [1/dt]
deltaepsl = 1.09


-- all important parameters were set above this point. If at a later stage there is need to change the properties of individual cylinders in an automated fashion, this can be accomplished by rewriting the following lines.

real_sy = (real_length-real_offset1-real_offsetN)/(n_cyl-1)


-- conversion to computation scale
conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

sx = 2*math.floor(real_sx/conv/2+.5)
sy = math.floor(real_sy/conv+.5)
sz = 2*math.floor(real_sz/conv/2+.5)
offset1 = math.floor(real_offset1/conv+.5)
offsetN = math.floor(real_offsetN/conv+.5)
rcyl = {}; hhcyl = {}; icyl = {}; jcyl = {}; kcyl = {};
for i = 1, n_cyl, 1 do
  rcyl[i] = math.floor(real_dcyl/2/conv+.5)
  hhcyl[i] = math.floor(real_hhcyl/conv+.5)
  icyl[i] = 0
  jcyl[i] = offset1+(i-1)*sy
  kcyl[i] = 0
end


imin = -sx/2
imax = sx/2
jmin = 0
jmax = offset1 + (n_cyl-1)*sy + offsetN
if (dim == 3) then 
  kmin = -sz/2
  kmax = sz/2
  dk = 3
elseif (dim == 2) then
  kmin = 0
  kmax = 0
  dk = 1
else
  printf("ERROR: Dimension not defined properly") 
end

jinj = math.floor(real_yinj/conv+.5)
jfft1 = math.floor(real_yfft1/conv+.5)
jfft2 = jmax - math.floor(real_yfft2/conv+.5)

-- print some parameters
print("Resolution:                                          ", resolution)
print("Conversion factor:                              dx = ", conv, "nm")
print("Courant factor:                                 dt = ", dt, "dx")
print("Wavelength (grid):                                   ", 1/inv_wavelength)
print("Inverse wavelength (grid):                           ", inv_wavelength)
print("Total length of structure: (grid)		    ", offset1+(n_cyl-1)*sy+offsetN) 
print("Offset between center of first cylinder and near PML (grid): ", offset1)
print("Offset between center of last cylinder and far PML (grid):   ", offset1+(n_cyl-1)*sy+offsetN-jcyl[n_cyl])
print("Distance between incoupling plane and first cylinder (grid): ", jcyl[1] - rcyl[1] - jinj)
print("Distance between reflection plane and first cylinder (grid): ", jcyl[1] - rcyl[1] - jfft1)
print("Distance between transmission plane and last cylinder (grid):", jfft2 - jcyl[n_cyl] - rcyl[n_cyl])
print("Distance of cylinders centers along k-vector (grid): ", sy)
print("Transverse distance of cylinder centers (grid):      ", sx)
print("Distance between cylinders (grid):                   ", sy - 2*math.floor(real_dcyl/2/conv+.5))
print("Number of cylinders:			            ", n_cyl)
for i,v in ipairs(rcyl) do
print("Cylinder ", i, ":")
print("Radius of cylinder (grid):             ", v)
if (dim == 3) then
  print("Half-Height of cylinder (grid):        ", hhcyl[i])
end
print("Position of center of cylinder (grid): ", icyl[i], jcyl[i], kcyl[i])
end

print("Computational window without PML:         ", imin, imax, jmin, jmax, kmin, kmax)
print("Time steps:                               ", ncycles)
max_tstep = math.floor(n_max*math.sqrt((imax-imin)^2+(jmax-jmin)^2+(kmax-kmin)^2)/dt +
            (attackl+sustainl+decayl)/inv_wavelength +.5)
print("Maximum number of time steps to")
print("propagate wave through structure:         ", max_tstep)
if (sy-2*math.floor(real_dcyl/2/conv+.5)<=1) then
   print("WARNING! Cylinders are very close and might overlap!")
end
