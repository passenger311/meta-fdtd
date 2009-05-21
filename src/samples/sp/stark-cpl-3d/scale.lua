
if TASKID then
   print("TASKID = ",TASKID)
   TASKPAR = 0.1 + (TASKID-1) * 0.05  -- vary distance between channels
else
   TASKPAR = 1.
end

--- select mode

te = true

if te then               -- te mode: Ey dominant

   betaeff = 2.5
   injfile = "tfsfey.set"
   injplane = { phi=0, theta=0, psi=0.0, nrefr=betaeff }

else                  -- tm mode: Ex dominant

   betaeff = 2.6
   injfile = "tfsfex.set"
   injplane = { phi=0, theta=0, psi=90.0, nrefr=betaeff }

end

--- geometrical parameters in real units(um)

real_wavelength   = 1.550               -- real wavelength 

real_length_wg1   = 2.4                 -- waveguide length
real_length_mmi   = 0.1                 -- mmi length
real_length_wg2   = 2.4                 -- waveguide length

real_width_mmi    = 10.0                -- mmi width 
real_width_wg     = 0.305               -- waveguide width   
real_hwidth_sep = 0.65                -- separation distance of two waveguide channels

real_height_wg    = 0.300               -- waveguide height
real_height_bsi   = 0.100               -- mmi height 
real_height_bsio2 = 2*real_height_wg    -- cladding height
real_height_top   = 0.100               -- upper cladding height

real_height = real_height_bsio2 +real_height_bsi + real_height_wg + real_height_top -- total height of structure
real_length = real_length_wg1 + real_length_mmi + real_length_wg2                   -- total length of the structure

angle = 10.           -- angle of rotation of the channels

--- material parameters

nref_si   = 3.47             -- refractive index of MMI and waveguide channels
nref_sio2 = 1.444            -- refractive index of the cladding

eps_si   = nref_si^2         -- dielectric constant of the MMI and waveguide channels
eps_sio2 = nref_sio2^2       -- dielectric constant of the base cladding
eps_bg = eps_sio2            -- dielectric constant of the background (top cladding)


--- scaling parameters

resolution = 30              -- resolution of wavelength in optically thickest medium (even number)         
dt         = 0.574           -- time step

real_dx       = real_wavelength / nref_si / resolution  -- conversion factor between real and grid length scale
invwavelength = real_dx / real_wavelength               -- this gives frequency for c=1, inverse of wavelength


--- waveguide effective refractive index (neff)

neff = 3.47


--- some other parameters
real_kfft0 = 0.1
real_kinj  = 0.2                               -- kinj distance (distance of injection plane from PML?)
real_kfft1 = real_length_wg1 - 0.4             -- fourier transform near to source
real_kfft3 = real_length     - 0.4             -- at the end of the channels before few grid points from PML
real_kfft2 = real_kfft3 - real_length_mmi      -- at the middle of the structure

--- in order to print do the following

print("resolution = ", resolution)
print("real_dx    = ", real_dx)
print("wavelength (grid) = ", real_wavelength/real_dx)
print("invwavelength (grid) = ", invwavelength)


--- Gaussian envelope of injection field 
--width1= 2   -- width in number of periods
--offset1 = 0 -- offset in number of periods
--attack1 = 4 -- attack in number of periods
--sustain1 = 0 -- sustain in number of periods
--decay1 = 4 -- decay in number of periods
--nrefr = 3.47 -- reference injection refractive index

--- time steps and excitation pulse

ncyc        = 6*1024       -- number of cycles 
pulsehwhm   = 150
pulsehsteps = 500          -- from center of the gaussian pulse to its tail on its both sides total is 1000


-- parameters in grid coordinates

height_wg    = math.floor(real_height_wg/real_dx   + 0.5)
height_bsi   = math.floor(real_height_bsi/real_dx  + 0.5)
height_bsio2 = math.floor(real_height_bsio2/real_dx+ 0.5)
height       = math.floor(real_height/real_dx      + 0.5)
hwidth_wg    = math.floor(real_width_wg/2/real_dx  + 0.5)
hwidth_mmi   = math.floor(real_width_mmi/2/real_dx + 0.5)
length_wg1   = math.floor(real_length_wg1/real_dx  + 0.5)
length_mmi   = math.floor(real_length_mmi/real_dx  + 0.5)
length       = math.floor(real_length/real_dx      + 0.5)
hwidth_sep   = math.floor(real_hwidth_sep/real_dx  + 0.5)
length_wg2   = math.floor(real_length_wg2/real_dx  + 0.5)
kinj         = math.floor(real_kinj/real_dx        + 0.5)
kfft0        = math.floor(real_kfft0/real_dx      + 0.5)
kfft1        = math.floor(real_kfft1/real_dx       + 0.5)
kfft2        = math.floor(real_kfft2/real_dx       + 0.5)
kfft3        = math.floor(real_kfft3/real_dx       + 0.5)
yc           = math.floor(height_bsio2+(height_bsi+height_wg)/2 + 0.5)

yc1          = math.floor(height_bsio2+ (height_bsi)/2+0.5)   -- for fourier transform lower point (y direction)
yc2          = math.floor(height_bsio2+ height_bsi+height_wg)   -- for fourier transform higher point  

--- print the following parameters

print("height_wg (grid)    = ", height_wg)
print("height_bsi (grid)   = ", height_bsi)
print("height_bsio2 (grid) = ", height_bsio2)
print("height (grid)       = ", height)
print("hwidth_wg (grid)    = ", hwidth_wg)
print("length_wg1 (grid)   = ", length_wg1)
print("length_mmi (grid)   = ", length_mmi)
print("hwidth_mmi (grid)   = ", hwidth_mmi)
print("hwidth_sep (grid)   = ", hwidth_sep)
print("length_wg2 (grid)   = ", length_wg2)
print("length (grid)       = ", length)
print("kinj (grid)         = ", kinj)
print("kfft0 (grid)        = ", kfft0)
print("kfft1 (grid)        = ", kfft1)
print("kfft2 (grid)        = ", kfft2)
print("yc (grid)           = ", yc)

print("yc1 (grid)           = ", yc1)
print("yc2 (grid)           = ", yc2)
--- pml cells

cpml = 11      -- in real units it will be real_dx times cpml 