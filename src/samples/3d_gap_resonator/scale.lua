if (taskid and (taskid ~= "undefined")) then
   print("taskid = ",taskid)
else
   taskid = 1 
end

pi = 3.141592653589793116


--- Real values in nm

real_wavelength = 650.0    -- real wavelength

real_hdist = 20/2           -- half distance of capacitor plates
real_hwidthx = 130/2      -- half width of capacitor plates
real_hwidthy = 130/2
real_hheight = 20/2          -- half height of capacitor plates (0 equals one yee cell)
real_redges = 4
real_hdist_tfsf_i = real_hwidthx + 4
real_hdist_tfsf_j = real_hwidthy + 4
real_hdist_tfsf_k = real_hdist + 2*real_hheight + 4


num_np = 0; -- number of nanoparticles

real_rnp = {}; real_xnp = {}; real_ynp = {}; real_znp = {}
math.randomseed( 1235 )

for i =1,num_np do
  real_rnp[i] = 2
  real_xnp[i] = math.random(-real_hwidthx + 1.5*real_rnp[i],real_hwidthx - 1.5*real_rnp[i])
  real_ynp[i] = math.random(-real_hwidthy + 1.5*real_rnp[i],real_hwidthy - 1.5*real_rnp[i])
  real_znp[i] = math.random(-real_hdist + 1.5*real_rnp[i],real_hdist - 1.5*real_rnp[i])
end


n_sphere = 3.
n_bg = 1.0
n_max = n_bg  -- maximum refractive index to determine optically thickest medium
mat = 'gold' -- gold, silver 

step_dft = 2
sampl_dft = 512
step_fft = 2
sampl_fft = 512

--- conversion to computation scale

resolution = 650/1 -- resolution of wavelength in optically thickest medium (even number)

conv = real_wavelength/resolution/n_max -- conversion factor between real and computation length scale
frequ_factor = 2.99792458e5  -- change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

inv_wavelength = conv/real_wavelength -- inverse wavelength

hdist = math.floor(real_hdist/conv+.5)
hwidthx = math.floor(real_hwidthx/conv+.5)
hwidthy = math.floor(real_hwidthy/conv+.5)
hheight = math.floor(real_hheight/conv+.5)
redges = math.floor(real_redges/conv+.5)
hdist_tfsf_i = math.floor(real_hdist_tfsf_i/conv+.5)
hdist_ntff_i = hdist_tfsf_i + 2
hdist_tfsf_j = math.floor(real_hdist_tfsf_j/conv+.5)
hdist_ntff_j = hdist_tfsf_j + 2
hdist_tfsf_k = math.floor(real_hdist_tfsf_k/conv+.5)
hdist_ntff_k = hdist_tfsf_j + 2


rnp={}; inp={}; jnp={}; knp={}
for i,v in ipairs(real_rnp) do
  rnp[i] = math.floor(real_rnp[i]/conv+.5)
  inp[i] = math.floor(real_xnp[i]/conv+.5)
  jnp[i] = math.floor(real_ynp[i]/conv+.5)
  knp[i] = math.floor(real_znp[i]/conv+.5)
end


--- Gaussian envelope of injection field
widthl = 1     -- width in number of periods
offsetl = 0    -- offset in number of periods
attackl = 4  -- attack in number of periods
sustainl = 0   -- sustain in number of periods
decayl = 8     -- decay in number of periods
nrefr = n_bg -- reference injection refractive index


--- PML size parameter

size_pml = 9


--- Padding size parameter

size_pad = 3


--- Courant factor

dt = 0.574  -- time step length compared to grid step length (--> Courant stability factor)

ncycles = 2*32768-1 -- number of cycles


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
print("ntff-domain x dir (grid):                 ", hdist_ntff_i)
print("tfsf-domain x dir (grid):                 ", hdist_tfsf_i)
print("ntff-domain y dir (grid):                 ", hdist_ntff_j)
print("tfsf-domain y dir (grid):                 ", hdist_tfsf_j)
print("ntff-domain z dir (grid):                 ", hdist_ntff_k)
print("tfsf-domain z dir (grid):                 ", hdist_tfsf_k)
print("Half distance of capacitor plates (grid): ", hdist)
print("Half width of capacitor plates, x (grid): ", hwidthx)
print("Half width of capacitor plates, y (grid): ", hwidthy)
print("Height of capacitor plates (grid):        ", 2*hheight+1)
print("Rounding of corners (grid):               ", redges)
print("Size of padding (grid):                   ", size_pad)
print("Size of PML (grid):                       ", size_pml)
print("The following parameters are only applicable in the case of nanoparticles")
for i,v in ipairs(rnp) do
  print("Radius of nanoparticle (grid):            ", v)
  print("Position of nanoparticle (grid):          ", inp[i], jnp[i], knp[i])
end


--- Write invlambda.in file for DFT

foutput = io.open("invlambda.in","w+")
foutput2 = io.open("lambda.in","w+")
for i = 600,1100,5 do
   foutput:write(conv/i, "\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("invlambda2.in","w+")
foutput2 = io.open("lambda2.in","w+")
for i = 600,1000,10 do
   foutput:write(conv/i,"\n")
   foutput2:write(i, "\n")
end
foutput:close()
foutput2:close()

foutput = io.open("data.save","w+")
foutput:write(conv,"\n")
foutput:write(2*(hdist_ntff_i+2), " ", 2*(hdist_ntff_j+2), " ", 2*(hdist_ntff_k+2), "\n")
foutput:write(2*(hdist_tfsf_i-2), " ", 2*(hdist_tfsf_j-2), " ", 2*(hdist_tfsf_k-2))
foutput:close()
