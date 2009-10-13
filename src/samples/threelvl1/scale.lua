
print("---- scale")


--- material parameters (um)

real_wavelength = 0.525   

eps_bg = 1.

nref = math.sqrt(eps_bg)

--- scaling parameters (in micrometers)

resolution = 30

dt = 0.99

real_dx       = real_wavelength / resolution 

invwavelength  = real_dx / real_wavelength


--- 

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)

