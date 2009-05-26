
if TASKID then
   print("TASKID = ",TASKID)
   TASKPAR = 0.1 + (TASKID-1) * 0.05  -- vary distance between channels
else
   TASKPAR = 1.
end

--- utility function; rounds to integer number

function round(f)
   return math.floor(f+0.5)
end


--- select mode

te = false

if te then               -- te mode: Ey dominant

   betaeff = 2.50
   injfile = "tfsfey.set"
   injplane = { phi=0, theta=0, psi=0, nrefr=betaeff }
   bc = 2 -- PMC

else                  -- tm mode: Ex dominant

   betaeff = 2.54
   injfile = "tfsfex.set"
   injplane = { phi=0, theta=0, psi=90, nrefr=betaeff }
   bc = 0 -- PEC

end

--- geometrical parameters in real units(um)

real_wavelength   = 1.550               -- real wavelength 

real_length_wg1   = 2.4                 -- waveguide length
real_length_mmi   = 0.1                 -- mmi length
real_length_wg2   = 2.4                 -- waveguide length
real_width_mmi    = 10.0                -- mmi width 
real_width_wg     = 0.305               -- waveguide width   
real_width_pad    = real_width_wg
real_hwidth_sep   = 0.65                -- separation (center to center)
real_height_wg    = 0.300               -- waveguide height
real_height_bsi   = 0.100               -- mmi height 
real_height_bsio2 = 2*real_height_wg    -- cladding height
real_height_top   = 0.100               -- upper cladding height

real_height = real_height_bsio2 +real_height_bsi + real_height_wg + real_height_top angle = 40.           -- angle of rotation of the channels

--- calculate alpha

alpha = math.pi/180. * angle

--- material parameters

nref_si   = 3.47             -- refractive index of MMI and waveguide channels
nref_sio2 = 1.444            -- refractive index of the cladding

eps_si   = nref_si^2         -- dielectric constant of the MMI and waveguide channels
eps_sio2 = nref_sio2^2       -- dielectric constant of the base cladding
eps_bg = eps_sio2            -- dielectric constant of the background (top cladding)


--- scaling parameters

resolution = 20              -- resolution of wavelength in optically thickest medium (even number)         
dt         = 0.574           -- time step

real_dx       = real_wavelength / nref_si / resolution  -- conversion factor between real and grid length scale
invwavelength = real_dx / real_wavelength               -- this gives frequency for c=1, inverse of wavelength

--- waveguide effective refractive index (neff)

neff = betaeff

--- print stuff

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)


-- parameters in grid coordinates

height_wg    = real_height_wg/real_dx
height_bsi   = real_height_bsi/real_dx  
height_bsio2 = real_height_bsio2/real_dx
height       = real_height/real_dx
hwidth_wg    = real_width_wg/2/real_dx  
hwidth_mmi   = real_width_mmi/2/real_dx 
width_pad    = real_width_pad/real_dx
length_wg1   = real_length_wg1/real_dx  
length_mmi   = real_length_mmi/real_dx  
hwidth_sep   = real_hwidth_sep/real_dx  
length_wg2   = real_length_wg2/real_dx  


--- calculate time steps 

total_length = length_mmi  + length_wg1  + length_wg2  
total_time = 2 * total_length * nref_si

ncyc        = round(total_time/dt)       -- number of time steps


--- excitation pulse

pulsehwhm   = 150
pulsehsteps = 500

--- calculate points for geometry!

l1 = length_wg1
l2 = length_mmi
l3 = length_wg2 * math.cos(alpha)
w1 = hwidth_mmi
w2 = hwidth_sep
w3 = hwidth_wg

alpha = math.pi/180. * angle

----

y1 = 0 - 100
y2 = l1
y3 = y2+l2
y4 = y3+l3 + 100

x1_1 = w3
x2_1 = x1_1
x2_2 = w1
x3_1 = w2
x3_2 = 2*w3/math.cos(alpha) + x3_1
x3_3 = w1
x4_1 = x3_1 + math.tan(alpha)*l3
x4_2 = x3_2 + math.tan(alpha)*l3

p1 = {x1_1,y1}
p2 = {x2_1,y2}
p3 = {x2_2 ,y2}
p4 = {x3_3,y3} 
p5 = {x3_2,y3}
p6 = {x4_2,y4}
p7 = {x4_1,y4}
p8 = {x3_1,y3}

p5 = { p5[1]-2*math.sin(alpha), p5[2]-2*math.cos(alpha) }
p8 = { p8[1]-2*math.sin(alpha), p8[2]-2*math.cos(alpha) }

points = {p1,p2,p3,p4,p5,p6,p7,p8}

--- calculate length and width of computational domain

hwidth = p6[1] + width_pad
length = p4[2] + length_wg2 * math.cos(alpha)


--- setup source / diagnostic planes for ffts

kinj = 10
kinf = kinj + 5
kinb = kinj - 5
kin0 = kinj + round(length_wg1/2)
kch  = round( length - 10 )
ich  = round( ( kch - y3 ) * math.tan(alpha) + hwidth_sep + hwidth_wg )

jc           = round(height_bsio2+(height_bsi+height_wg)/2)
jc1          = round(height_bsio2+ (height_bsi)/2 )
jc2          = round( height_bsio2+ height_bsi+height_wg)
jstep = 2

idelta = round(2*hwidth_wg)
istep = 2

--- print the following parameters

print("height_wg    = ", height_wg)
print("height_bsi   = ", height_bsi)
print("height_bsio2 = ", height_bsio2)
print("height       = ", height)
print("hwidth_wg    = ", hwidth_wg)
print("hwidth_mmi   = ", hwidth_mmi)
print("hwidth_sep   = ", hwidth_sep)
print("hwidth       = ", hwidth)
print("length_wg1   = ", length_wg1)
print("length_mmi   = ", length_mmi)
print("length_wg2   = ", length_wg2)
print("length       = ", length)
print("kinj         = ", kinj)
print("kinf         = ", kinf)
print("kinb         = ", kinb)
print("kin0         = ", kin0)
print("kch          = ", kch)
print("ich          = ", ich)

print("jc           = ", jc)
print("jc1          = ", jc1)
print("jc2          = ", jc2)

for i,v in ipairs(points) do
   print("P"..tostring(i).." = ",v[1],v[2])
end

--- pml cells

cpml = 11      -- in real units it will be real_dx times cpml 
