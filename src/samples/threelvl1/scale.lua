
print("---- scale")


--- material parameters (um)

real_wavelength = 0.525   

epsbg = 4.

nref = math.sqrt(epsbg)

--- scaling parameters (in micrometers)

resolution = 50

dt = 0.999

real_dx       = real_wavelength / resolution 

invwavelength  = 1 / resolution

pi = 3.14159265

mu = 0.01

tau = 1000
area = 1 --- in units of 2 pi
unitCharge = 1.602176487e-19

hbar = 1.05457148e-34
clight = 299792458
eps0 = 8.85418781e-12

tauSI = tau * real_dx * 1.e-6 * dt / clight
muSI = mu * real_dx * 1.e-6 * unitCharge
ESI = 2 * area * hbar / muSI / tauSI 
ENHL = ESI * real_dx * math.sqrt(real_dx) * 1.e-9 * math.sqrt(eps0) / clight * math.log(2+math.sqrt(3))


--- 

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)

