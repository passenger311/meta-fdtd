
print("---- scale")

--- material parameters

nref = {}

nref.air  = 1
nref.si   = 3.2 
nref.sio2 = 1.444
nref.bg = nref.sio2

eps = {}

for k,v in pairs(nref) do 
   eps[k] = v^2
end


--- scaling parameters

real_wavelength = 1.55
resolution = 50     
dt         = 0.705

real_dx       = real_wavelength / nref.si / resolution 
invwavelength = real_dx / real_wavelength 



--- in order to print do the following

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)

