
taskid = select(1,...)

if taskid then
   print("taskid = ",taskid)
-- do something like this:
   scanval = 4. + (taskid-1) * 0.5
-- for array job set real_length_mmi=scanval.
end


--- real device parameters (um)

real_wavelength = 1.550
real_width_mmi = 2.43
real_length_mmi = 6.15
real_width_wg = 0.400
real_height_wg = 0.220
real_length_wg1 = 2.4
real_length_wg2 = 2.4
real_hwidth_sep = 0.65
real_kinj = 0.2
real_height_bsio2 = 2*real_height_wg
real_height = real_height_bsio2 + real_height_wg + real_height_wg
real_length = real_length_wg1 + real_length_mmi + real_length_wg2


real_kfft1 = real_kinj + 1.6 
real_kfft3 = real_length - 0.4
real_kfft2 = real_kfft3 - real_length_mmi

--- material parameters

nref_si = 3.47
nref_sio2 = 1.444

eps_si = nref_si^2
eps_sio2 = nref_sio2^2

--- scaling parameters

resolution = 45
dt = 0.574

real_dx = real_wavelength / nref_si / resolution
invwavelength = real_dx / real_wavelength

print("resolution = ", resolution)
print("real_dx = ", real_dx)
print("wavelength (grid) = ", real_wavelength/real_dx)
print("invwavelength (grid) = ", invwavelength)



--- time steps and excitation pulse

ncyc = 8192+4096
pulsehwhm = 150
pulsehsteps = 500


-- parameters in grid coordinates

height_wg = math.floor(real_height_wg/real_dx+0.5)
height_bsio2 = math.floor(real_height_bsio2/real_dx+0.5)
height = math.floor(real_height/real_dx+0.5)
hwidth_wg = math.floor(real_width_wg/2/real_dx+0.5)
hwidth_mmi = math.floor(real_width_mmi/2/real_dx+0.5)
length_wg1 = math.floor(real_length_wg1/real_dx+0.5)
length_mmi = math.floor(real_length_mmi/real_dx+0.5)
length = math.floor(real_length/real_dx+0.5)
hwidth_sep = math.floor(real_hwidth_sep/real_dx+0.5)
length_wg2 = math.floor(real_length_wg2/real_dx+0.5)
kinj =  math.floor(real_kinj/real_dx+0.5)
kfft1 =  math.floor(real_kfft1/real_dx+0.5)
kfft2 =  math.floor(real_kfft2/real_dx+0.5)
kfft3 =  math.floor(real_kfft3/real_dx+0.5)
yc = height_bsio2+math.floor(height_wg/2+0.5)

print("height_wg (grid) = ", height_wg)
print("height_bsio2 (grid) = ", height_bsio2)
print("height (grid) = ", height)
print("hwidth_wg (grid) = ", hwidth_wg)
print("length_wg1 (grid) = ", length_wg1)
print("length_mmi (grid) = ", length_mmi)
print("hwidth_mmi (grid) = ", hwidth_mmi)
print("hwidth_sep (grid) = ", hwidth_sep)
print("length_wg2 (grid) = ", length_wg2)
print("length (grid) = ", length)
print("kinj (grid) = ", kinj)
print("kfft1 (grid) = ", kfft1)
print("kfft2 (grid) = ", kfft2)
print("yc (grid) = ", yc)


--- pml cells

cpml = 11