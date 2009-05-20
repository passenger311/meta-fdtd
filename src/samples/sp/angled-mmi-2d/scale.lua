
--- This is like the 3d version but with k<->j flipped

if TASKID then
   print("TASKID = ",TASKID)
   TASKPAR = 1.0 + (TASKID-1) * 0.05
else
   TASKPAR = 1.0                     
end

--- utility function; rounds to integer number

function round(f)
   return math.floor(f+0.5)
end


--- select mode

te = true

if te then               -- te mode: Ey dominant

   betaeff = 2.532
   injfile = "tfsfey.set"
   injplane = { phi=90, theta=90, psi=90.0, nrefr=betaeff }

else                  -- tm mode: Ex dominant

   betaeff = 2.530
   injfile = "tfsfex.set"
   injplane = { phi=90, theta=90, psi=0.0, nrefr=betaeff }

end

--- geometrical parameters in real units (um)

numchannel = 3

real_wavelength   = 1.550                -- real wavelength 

real_hlength_mmi  = 12.                 -- mmi half length
real_hwidth_mmi   = 1.0                 -- mmi width 
real_length_wg    = 1.5                 -- mmi length

real_width_wg     = 0.305               -- waveguide width   
real_width_pad    = real_width_wg/3. 
real_width_sep    = 0.305               -- separation distance of two waveguide channels

real_height_wg    = 0.300               -- waveguide height
real_height_bsi   = 0.100               -- mmi height 
real_height_bsio2 = 1.0*real_height_wg  -- cladding height
real_height_top   = 0.200               -- upper cladding height

real_height = real_height_bsio2 +real_height_bsi + real_height_wg + real_height_top -- total height of structure

angle = 10.           -- angle of rotation of mmi axis vs input channel (alpha)

--- material parameters

nref_si   = 3.47             -- refractive index of MMI and waveguide channels
nref_sio2 = 1.444            -- refractive index of the cladding

eps_si    = nref_si^2        -- dielectric constant of the MMI and waveguide channels
eps_sio2  = nref_sio2^2      -- dielectric constant of the base cladding
eps_bg    = eps_sio2         -- dielectric constant of the background (top cladding)


--- scaling parameters

resolution = 20              -- resolution of wavelength in optically thickest medium (even number)         
dt         = 0.7             -- time step

real_dx       = real_wavelength / nref_si / resolution  -- conversion factor between real and grid length scale
invwavelength = real_dx / real_wavelength               -- this gives frequency for c=1, inverse of wavelength


--- waveguide effective refractive index (neff)

neff = betaeff


--- in order to print do the following

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)

--- time steps and excitation pulse

ncyc        = 6*1024       -- number of cycles 
pulsehwhm   = 150
pulsehsteps = 500          -- from center of the gaussian pulse to its tail on its both sides total is 1000


-- parameters in grid coordinates

height_wg    = real_height_wg/real_dx
height_bsi   = real_height_bsi/real_dx
height_bsio2 = real_height_bsio2/real_dx
height       = real_height/real_dx
hwidth_wg    = real_width_wg/2/real_dx
hwidth_mmi   = real_hwidth_mmi/real_dx
width_pad    = real_width_pad/real_dx
width_sep    = real_width_sep/real_dx
hlength_mmi  = real_hlength_mmi/real_dx
length_wg    = real_length_wg/real_dx

yc           = height_bsio2+ (height_bsi+height_wg)/2
yc1          = height_bsio2+ (height_bsi)/2   -- for fourier transform lower point (y direction)
yc2          = height_bsio2+ height_bsi+height_wg   -- for fourier transform higher point  


--- precalculate

l = 100000000 -- (something really long!)
l1 = hlength_mmi
w2 = hwidth_wg
w1 = hwidth_mmi

ys = height_bsio2 + height_bsi
ye = ys + height_wg

alpha = math.pi/180. * angle
beta = math.pi/180. * ( 90.-angle )

--- check whether all out channels fit

real_mmi_xproj = math.sin(alpha)*real_hlength_mmi*2
real_ch_xproj = numchannel*real_width_wg + (numchannel-1)*real_width_sep

if ( real_ch_xproj > real_mmi_xproj ) then
   error("output channels won't fit onto mmi")
end

--- calculate points of angled mmi

x1 = l1*math.sin(alpha)
z1 = l1*math.cos(alpha)
x2 = -x1
z2 = -z1
x3 = x1 - w1*math.sin(beta)
z3 = z1 + w1*math.cos(beta)
x4 = -x3
z4 = -z3
x5 = x1 + w1*math.sin(beta)
z5 = z1 - w1*math.cos(beta)
x6 = -x5
z6 = -z5
x7 = x4
z7 = z4 + 2*w2 / math.tan(beta)
x8 = x4 - 2*w2
z8 = - l
x9 = -x7
z9 = z3 - 2*w2 / math.tan(alpha)
x10 = x3 - 2*w2
z10 = l

ch_dx = width_sep + 2*w2
ch_dz = ch_dx / math.tan(alpha)

p1 = { x1,z1 }
p2 = { x2,z2 }
p3 = { x3,z3 }
p4 = { x4,z4 }
p5 = { x5,z5 }
p6 = { x6,z6 }
p7 = { x7,z7 }
p8 = { x8,z8 }
p9 = { x9,z9 }
p10 = { x10,z10 }

points = {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10}

--- calculate length and width of computational domain

hwidth = math.max(x5,x6) + width_pad
hlength = z3 + length_wg

--- setup source / diagnostic planes

jinj =  round(- hlength + 20)
jfft1 = round(- hlength + 10)
jfft2 = round(- hlength + length_wg / 2)
jfft3 = round(hlength - 10)

--- x pos of channel centers

iin =  round(x8 + hwidth_wg)
iout1 = round(x10 + hwidth_wg)
iout2 = round(iout1 - width_sep - 2*hwidth_wg)
iout3 = round(iout2 - width_sep - 2*hwidth_wg)
iout4 = round(iout3 - width_sep - 2*hwidth_wg)
iout5 = round(iout4 - width_sep - 2*hwidth_wg)
iout6 = round(iout5 - width_sep - 2*hwidth_wg)

deltai = round(hwidth_wg+width_sep/2)

kc = round(yc)

kc1 = round(yc1)
kc2 = round(yc2)

--- print the following parameters

print("height_wg     = ", height_wg)
print("height_bsi    = ", height_bsi)
print("height_bsio2  = ", height_bsio2)
print("height        = ", height)
print("hwidth        = ", hwidth)
print("hwidth_wg     = ", hwidth_wg)
print("hwidth_mmi    = ", hwidth_mmi)
print("width_sep     = ", width_sep)
print("hlength       = ", hlength)
print("length_wg     = ", length_wg)
print("hlength_mmi   = ", hlength_mmi)
print("kc            = ", kc)
print("kc1           = ", kc1)
print("kc2           = ", kc2)

print("jinj          = ",jinj)
print("jfft1         = ",jfft1)
print("jfft2         = ",jfft2)
print("jfft3         = ",jfft3)

for i,v in ipairs(points) do
   print("P"..tostring(i).." = ",v[1],v[2])
end

--- pml cells

cpml = 11     -- in real units it will be real_dx times cpml 

