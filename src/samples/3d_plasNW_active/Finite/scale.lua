pi = 3.141592653589793116

--- Real values in nm

real_cyl_radius = taskid/2--55.0 --radius of silver layers
real_hspacer = 10.0 --radius of spacer layer
real_cyl_length = 1600.0 --y dimension of perforation
real_hmetal_sub = 130.0 --z dimension radius of metal substrate
real_dist = 600 --distance of measuring plane from structure
real_wavelength = 300.0    -- real wavelength

n_bg = 1.0
n_spacer = 1.4
n_substrate = 1.0
n_nanowire = 2.3
mat = 'silver' -- gold, silver, dielectric with n_perf

step_fft = 10  -- spatial step of fft
step_dft = 1  -- spatial step of dft
samp_fft = 8192*8/4 -- sampling of fft
samp_dft = 8192*8/4 --sampling of dft

--- conversion to computation scale

resolution = real_wavelength/10 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

cyl_radius = real_cyl_radius/conv --math.floor(real_cyl_radius/conv+.5)
cyl_diameter = 2*cyl_radius
hspacer = math.floor(real_hspacer/conv+.5)
cyl_length = math.floor(real_cyl_length /2 /conv+.5)
hmetal_sub = math.floor(real_hmetal_sub/conv+.5)
dist = math.floor(real_dist/conv+.5)


--- PML size parameter
size_pml = 12

--- Padding size parameter
size_pad = 3


hdist_tfsf_i = math.ceil(cyl_radius)+dist --2*cyl_diameter 
hdist_ntff_i = hdist_tfsf_i
hdist_tfsf_j = cyl_length+dist --2*cyl_diameter
hdist_ntff_j = hdist_tfsf_j 
hdist_tfsf_kn = hmetal_sub
hdist_ntff_kn = hdist_tfsf_kn
hdist_tfsf_kp = cyl_diameter + hspacer + dist
hdist_ntff_kp = hdist_tfsf_kp

--- Gaussian envelope of injection field
pump_on = false
inv_wl = inv_wavelength
ampl = 5e6   -- electric field strength in V/m
widthl = 400    -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 0--1200  -- attack in number of periods
sustainl = 0--200*(taskid-1)  -- sustain in number of periods
decayl = 0--4000     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index
psi = 90

probe_on = false
pinv_wl = conv/380
pampl = 5e5 --1e7   -- electric field strength in V/m
pwidthl = 1    -- width in number of periods
poffsetl = 0    -- offset in number of periods
pattackl = attackl+sustainl+decayl+10 -- attack in number of periods
psustainl = 0   -- sustain in number of periods
pdecayl = 200    -- decay in number of periods
ppsi = 0

--- Fourlevel system parameters
invlambdala = pinv_wl			-- 380 nm
invlambdalb = inv_wl			-- 300 nm
pump_rate = 1*conv/frequ_factor		-- 1 per ps
gammala = 1/8*1000*conv/frequ_factor	-- 4 fs
gammalb = 1/8*1000*conv/frequ_factor	-- 4 fs
dipolea = .2/conv			-- 0.1 nm
dipoleb = .2/conv			-- 0.1 nm
density = 0.1*conv^3			-- 6e18 1/cm^3
init_inv = 0.7
gamma30 = 0
gamma32 = 1000/100*conv/frequ_factor	-- 100 fs
gamma21 = 1/500*conv/frequ_factor	-- 500 ps
gamma10 = 1000/100*conv/frequ_factor	-- 100 fs
localfield = 0
volfactor = 10000
linefactor = 1

--- diagnostics?
diagebal=true
diagpspec=true
diagmode=true


--- Courant factor
dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles_probe = 8*16384-1
ncyc_probe_start = 2*16384--math.floor((attackl+sustainl+decayl)*resolution/dt+.5)
ncycles = ncyc_probe_start+ncycles_probe -- number of cycles


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
  real_gammaDL = 12.6126e-12
  real_omegaL = 8257.18/2/pi
  real_gammaL = 195.614 --195.614
  deltaepsl = 2.23994
  real_omegaL2 = 3057.07/2/pi 
  real_gammaL2 = 852.675 --852.675
  deltaepsl2 = 0.222651
elseif (mat == 'dielectric') then
  eps_infDL = n_perf^2
end

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("ntff-domain x dir (grid):                 ", hdist_ntff_i)
print("tfsf-domain x dir(grid):                  ", hdist_tfsf_i)
print("ntff-domain y dir (grid):                 ", hdist_ntff_j)
print("tfsf-domain y dir(grid):                  ", hdist_tfsf_j)
print("ntff-domain pos z dir (grid):             ", hdist_ntff_kp)
print("ntff-domain neg z dir (grid):             ", hdist_ntff_kn)
print("tfsf-domain pos z dir (grid):             ", hdist_tfsf_kp)
print("tfsf-domain neg z dir (grid):             ", hdist_tfsf_kn)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
print("radius of metal (grid):                   ", cyl_radius)
print("radius of spacer (grid):                  ", hspacer)
print("Half-length of perforation in ydir (grid):", cyl_length)
print("distance of tfsf-domain from struct(grid):", dist)



--- Write invlambda.in file for DFT

foutput = io.open("invlambda2.in","w+")
foutput2 = io.open("lambda2.in","w+")
--[[
for i = 375,385,1 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
--]]

i=380
foutput:write(conv/i,"\n")
foutput2:write(i, "\n")
--[[
i=876
foutput:write(conv/i,"\n")
foutput2:write(i, "\n")
--]]
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(hspacer+math.ceil(cyl_radius),"\n")
foutput:write(n_bg,"\n")
foutput:write(n_substrate,"\n")
foutput:write(inv_wavelength,"\n")
foutput:write(hdist_tfsf_kn,"\n")
foutput:write(hdist_tfsf_kp,"\n")

foutput:close()
