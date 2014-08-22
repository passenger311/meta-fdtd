pi = 3.141592653589793116

--- Real values in nm

real_cyl_radius = cyl_rad	-- radius of ZnO wire
real_hspacer = 10.0		-- thickness of spacer layer
real_cyl_length = 400.0		-- length of cylinder
real_hmetal_sub = 100.0		-- thickness of substrate
real_dist = distance--2000			-- distance of top, sides end end of wire
real_disth = 450
real_wavelength = lam--390		-- real wavelength

n_refr = 1.8

n_bg = 1.0
n_nanowire = 2.3
mat = 'silver' -- gold, silver, silver2 (Drude only), dielectric with n_diel
n_diel = 1.0
if (mat=='dielectric') then
  n_spacer = 1.0
else
  n_spacer = 1.4
end

-- reflection setup?
if (refl=='1') then
  print('Reflection setup!')
  mult_factor=10 -- reflection setup for mult_factor = 1; mult_factor = 10 -> infinite wire
else
  mult_factor=100
end

step_fft = 1  -- spatial step of fft
step_dft = 1  -- spatial step of dft
samp_fft = 8192 -- sampling of fft
samp_dft = 8192*4 --sampling of dft

--- conversion to computation scale

resolution = real_wavelength/2 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

cyl_radius = math.floor(real_cyl_radius/conv+.5)
cyl_radius1 = real_cyl_radius/conv

cyl_diameter = 2*cyl_radius
hspacer = math.floor(real_hspacer/conv+.5)
cyl_length = math.floor(real_cyl_length /2 /conv+.5)
hmetal_sub = math.floor(real_hmetal_sub/conv+.5)
dist = math.floor(real_dist/conv+.5)
disth = math.floor(real_disth/conv+.5)


--- PML size parameter
size_pml = 10

--- Padding size parameter
size_pad = 3


hdist_tfsf_i = dist 
hdist_tfsf_jn = hmetal_sub
hdist_tfsf_jp = hspacer + disth
hdist_tfsf_kn = cyl_length+dist 
hdist_tfsf_kp = cyl_length+dist


--- Geometry processing for meta?
geo_on = false
prev_on = false

--- VTK output?
field_vtk_on = true
injection_vtk_on = false

--- Gaussian envelope of injection field
pump_on = false
inv_wl = inv_wavelength
ampl = 6e6   -- electric field strength in V/m
widthl = 100    -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 0--500--1200  -- attack in number of periods
sustainl = 0--200*(taskid-1)  -- sustain in number of periods
decayl = 0--4000     -- decay in number of periods
nrefr = n_refr -- reference injection refractive index
psi = 90

probe_on = false
probe_on_injx = true
probe_on_injy = true
pinv_wl = conv/lam--390
pampl = 1.0075029765e26 --1e7   -- electric field strength in V/m
pwidthl = 40--0    -- width in number of periods
poffsetl = 0    -- offset in number of periods
pattackl = (attackl+sustainl+decayl)*pinv_wl/inv_wl+160--300 -- attack in number of periods
psustainl = 20000   -- sustain in number of periods
pdecayl = 20000    -- decay in number of periods
ppsi = 90

--- Fourlevel system parameters
on_4lvl = false
invlambdala = pinv_wl			-- 380.1 nm
invlambdalb = inv_wl			-- 300 nm
pump_rate = 1*conv/frequ_factor		-- 1 per ps
gammala = 1/8*1000*conv/frequ_factor	-- 8 fs
gammalb = 1/8*1000*conv/frequ_factor	-- 8 fs
dipolea = .2/conv			-- 0.2 nm
dipoleb = .2/conv			-- 0.2 nm
density = 0.04*conv^3			-- 4e19 1/cm^3
init_inv = .5
gamma30 = 0
gamma32 = 1000/100*conv/frequ_factor	-- 100 fs
gamma21 = 1/500*conv/frequ_factor	-- 500 ps
gamma10 = 1000/100*conv/frequ_factor	-- 100 fs
localfield = 0
volfactor = 1000
linefactor = 1

--- diagnostics?
diagebal=true
diagpspec=false
diagmode=true


--- Courant factor
dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles_probe = 4*16380.1-1
ncyc_probe_start=math.floor((2*pattackl)*resolution/dt+.5)
--math.floor((pattackl+1.5*pwidthl)/inv_wl/dt+0.5)--0--math.floor((2*pattackl)*resolution/dt+.5)
--ncyc_probe_start = 10000--8*16380.1 --math.floor((attackl+sustainl+decayl)*resolution/dt+.5)
ncycles = ncyc_probe_start+ncycles_probe -- number of cycles


--- Drude-Lorentz material in THz

if (mat == 'gold') then
  eps_infDL = 5.9673
  real_omegaDL = 2113.6      -- Drude plasma frequency [2 pi c]
  real_gammaDL = 15.92*2*pi   -- Drude damping constant [1/dt]
  real_omegaL = 650.07       -- Lorentzian plasma frequency [2 pi c]
  real_gammaL = 104.86/2*2*pi   -- Lorentzian damping constant [1/dt]
  deltaepsl = 1.09
elseif (mat == 'silver2') then
  eps_infDL = 6.54--4.494--1.05
  real_omegaDL = 15600/2/pi--13910/2/pi -- plasma frequency in THz 
  real_gammaDL = 0--4*32.3 -- loss rate in metal in E12 rads/sec
elseif (mat == 'silver') then
  eps_infDL = 1.17152
  real_omegaDL = 13960.4/2/pi
  real_gammaDL = 12.6126e-12
  real_omegaL = 8257.18/2/pi
  real_gammaL = 195.614--/100
  deltaepsl = 2.23994
  real_omegaL2 = 3057.07/2/pi 
  real_gammaL2 = 852.675--/100
  deltaepsl2 = 0.222651
elseif (mat == 'dielectric') then
  eps_infDL = n_diel^2
end

--- print some parameters
print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("tfsf-domain x dir(grid):                  ", hdist_tfsf_i)
print("tfsf-domain pos y dir(grid)               ", hdist_tfsf_jp)
print("tfsf-domain neg y dir(grid)               ", hdist_tfsf_jn)
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
for i = 360,400,2.5 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
--]]
foutput:write(pinv_wl,"\n")
foutput2:write(conv/pinv_wl, "\n")
--[[
i=876
foutput:write(conv/i,"\n")
foutput2:write(i, "\n")
--]]
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(hspacer+cyl_radius,"\n")
foutput:write(n_bg,"\n")
foutput:write(inv_wavelength,"\n")
foutput:write(hdist_tfsf_kn,"\n")
foutput:write(hdist_tfsf_kp,"\n")

foutput:close()
