if (taskid and (taskid ~= "undefined")) then
   print("taskid = ",taskid)
else
   taskid = 1 
end

pi = 3.141592653589793116


--- Real values in nm

real_wavelength = 650.0    -- real wavelength

real_hdist = 20/2          -- half distance of plates
real_hwidth = 130/2        -- half width of plates
real_hheight = 20/2        -- half height of plates
real_redges = 4             -- radius of rounded edges
real_hdist_tfsf_i = real_hwidth + 20
real_hdist_tfsf_j = 2*real_hheight + real_hdist + 20
real_hdist_ntff_i = real_hdist_tfsf_i + 22
real_hdist_ntff_j = real_hdist_tfsf_j + 22 

n_sphere = 3.
n_bg = 1.0
n_max = n_bg  -- maximum refraxctive index to determine optically thickest medium
mat = 'silver' -- gold, silver 

step_dft = 1
sampl_dft = 512
step_fft = 1
sampl_fft = 512


--- conversion to computation scale

resolution = 650/.5 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hdist = math.floor(real_hdist/conv+.5)
hwidth = math.floor(real_hwidth/conv+.5)
hheight = math.floor(real_hheight/conv+.5)
redges = math.floor(real_redges/conv+.5)
hdist_tfsf_i = math.floor(real_hdist_tfsf_i/conv+.5)
hdist_ntff_i = math.floor(real_hdist_ntff_i/conv+.5)
hdist_tfsf_j = math.floor(real_hdist_tfsf_j/conv+.5)
hdist_ntff_j = math.floor(real_hdist_ntff_j/conv+.5)


--- Gaussian envelope of injection field

widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index


--- PML size parameter

size_pml = 12


--- Padding size parameter

size_pad = 3


--- Courant factor

dt = .7 -- 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 4*32768-1 -- number of cycles


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

eps_C_infDL = 2.554
real_C_omegaDL = 2149.77/2/pi
real_C_gammaDL = 9475.68
real_C_omegaL = 6820.0/2/pi
real_C_gammaL = 4799.37
C_deltaepsl = 3.56779
real_C_omegaL2 = 2833.44/2/pi 
real_C_gammaL2 = 1826.16
C_deltaepsl2 = 2.51021

eps_diel = n_sphere^2


--- print some parameters

print("Resolution:                               ", resolution)
print("Conversion factor:                   dx = ", conv, "nm")
print("Courant factor:                      dt = ", dt, "dx")
print("Wavelength (grid):                        ", 1/inv_wavelength)
print("Inverse wavelength (grid):                ", inv_wavelength)
print("ntff-domain x/y dir (grid):               ", hdist_ntff_i)
print("tfsf-domain x/y dir(grid):                ", hdist_tfsf_i)
print("ntff-domain z dir (grid):                 ", hdist_ntff_j)
print("tfsf-domain z dir (grid):                 ", hdist_tfsf_j)
print("Half distance of strips (grid):           ", hdist)
print("Half width of strips (grid):              ", hwidth)
print("Height of strips (grid):                  ", 2*hheight+1)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)


--- Write output files for DFT

foutput = io.open("invlambda.in","w+")
foutput2 = io.open("lambda.in","w+")
for i = 500,900,2 do
   foutput:write(conv/i, "\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("invlambda2.in","w+")
foutput2 = io.open("lambda2.in","w+")
for i = 600,800,5 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(2*(hdist_ntff_i+2), " ", 2*(hdist_ntff_j+2), " ", 1, "\n")
foutput:write(2*(hdist_tfsf_i-2), " ", 2*(hdist_tfsf_j-2), " ", 1)
foutput:close()
