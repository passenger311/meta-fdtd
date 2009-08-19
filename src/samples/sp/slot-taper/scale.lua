
--[[if TASKID then
   print("TASKID = ",TASKID)
   TASKPAR = 0.1 + (TASKID-1) * 0.05  -- vary distance between channels
else
   TASKPAR = 1.
end--]]

--- utility function; rounds to integer number

function round(f)
   return math.floor(f+0.5)
end


--- select mode

tm = false

if tm then               -- te mode: Ey dominant

   betaeff = 2.50
   injfile = "tfsfey.set"
   injplane = { phi=0, theta=0, psi=0, nrefr=betaeff }
   bc = 2 -- PMC

else                  -- te mode: Ex dominant

   betaeff = 2.7
   injfile = "tfsfex.set"
   injplane = { phi=0, theta=0, psi=90, nrefr=betaeff }
   bc = 0 -- PEC

end

--- geometrical parameters in real units(um)

cladding = true                          -- false=nometalcover, true=metal cover

real_wavelength   = 1.550                -- real wavelength 

real_width_edw    = 0.2                  -- evanescent decay width

real_length_i1     = 2.0                 -- length in 1
real_length_i2     = 0.5                 -- 0.5length in 2
real_length_it     = 5.0                 -- 2.5length in taper
real_length_c      = 1.0                 -- 2.0length channel 
real_length_ot     = 2.0                 -- 2.0length out taper
real_length_o2     = 0.5                 -- 0.5length out 2
real_length_o1     = 2.0                 -- 1.0length out 1
real_width_i       = 3.0                 -- width in 
real_width_c       = 0.0                 -- 0.400width channel
real_width_o       = 0.4                 -- width out 
real_width_ii      = 0.10               -- width clad in 
real_width_cc      = 0.260               -- width clad channel
real_width_oo      = 0.10                -- 0.2width clad out
real_width_gi        = 0.05                -- 0.05width sio2 seperation gap
real_width_go       = 0.05                -- 0.05width sio2 seperation gap
real_width_gc       = 0.025                -- 0.05width sio2 seperation gap
real_height_g      = 0.05                -- height sio2 seperation gap
real_width_pad     = 0.5

real_height_wg     = 0.200               -- waveguide height
real_height_bsi    = 0.000               
real_height_pad    =  2*real_height_wg   -- upper/lower cladding height


--[[
if not cladding then

real_width_ii      = 0.
real_width_cc      = 0.
real_width_oo      = 0.
real_width_gi        = 0.
real_width_go       = 0.
real_width_gc       = 0.

end
]]--

real_length_i = real_length_i1 + real_length_i2
real_length_o = real_length_o1 + real_length_o2
real_length_t = real_length_it + real_length_c + real_length_ot

real_height = real_height_pad +real_height_bsi + real_height_wg + real_height_pad
real_width = math.max(real_width_i+real_width_gi , real_width_c + real_width_gc, real_width_o + real_width_go) + real_width_pad
real_length = real_length_i + real_length_t + real_length_o

--- angles

alpha = math.atan(  ( real_width_i/2. + real_width_gi - real_width_c/2. - real_width_gc ) / real_length_it )
alpha2 = math.atan(  ( real_width_i/2. + real_width_ii + real_width_gi - real_width_c/2. - real_width_cc - real_width_gc ) / real_length_it )

beta = math.atan(   ( real_width_o/2. + real_width_go - real_width_c/2.  - real_width_gc ) /real_length_ot )
beta2 = math.atan(   ( real_width_o/2. + real_width_oo + real_width_go - real_width_c/2. -real_width_cc  - real_width_gc ) /real_length_ot )


--- material parameters

nref_si   = 3.47             -- refractive index of MMI and waveguide channels
nref_sio2 = 1.444            -- refractive index of the cladding

eps_si   = nref_si^2         -- dielectric constant of the MMI and waveguide channels
eps_sio2 = nref_sio2^2       -- dielectric constant of the base cladding
eps_bg = eps_sio2            -- dielectric constant of the background (top cladding)


--- scaling parameters

resolution = 35              -- resolution of wavelength in optically thickest medium (even number)         
dt         = 0.574           -- time step

real_dx       = real_wavelength / nref_si / resolution  -- conversion factor between real and grid length scale
invwavelength = real_dx / real_wavelength               -- this gives frequency for c=1, inverse of wavelength

--- waveguide effective refractive index (neff)
nrefr = 3.47
neff = betaeff 

--- print stuff

print("resolution    = ", resolution)
print("real_dx       = ", real_dx)
print("wavelength    = ", real_wavelength/real_dx)
print("invwavelength = ", invwavelength)


-- parameters in grid coordinates

length_i1    = real_length_i1/real_dx        -- length in 1
length_i2    = real_length_i2/real_dx        -- length in 2
length_it    = real_length_it/real_dx        -- length in taper
length_c     = real_length_c/real_dx         -- length channel 
length_ot    = real_length_ot/real_dx        -- length out taper
length_o2    = real_length_o2/real_dx        -- length out 2
length_o1    = real_length_o1/real_dx        -- length out 1
hwidth_i     = real_width_i/real_dx/2.       -- width in 
hwidth_c     = real_width_c/real_dx/2.       -- width channel
hwidth_o     = real_width_o/real_dx/2.       -- width out
width_ii     = real_width_ii/real_dx         -- width in 
width_cc     = real_width_cc/real_dx         -- width channel
width_oo     = real_width_oo/real_dx         -- width out
width_gi      = real_width_gi/real_dx          -- width sio2 seperation gap
width_go      = real_width_go/real_dx          -- width sio2 seperation gap
width_gc     = real_width_gc/real_dx          -- width sio2 seperation gap
height_g     = real_height_g/real_dx         -- height sio2 seperation gap
width_pad    = real_width_pad/real_dx
width_edw    = real_width_edw/real_dx


height_wg    = real_height_wg/real_dx        -- waveguide height
height_bsi   = real_height_bsi/real_dx       -- mmi height 
height_pad   = real_height_pad/real_dx       -- upper/lower cladding height

height_bsio2 = height_pad

length_i = real_length_i/real_dx
length_o = real_length_o/real_dx
length_t = real_length_t/real_dx

height = real_height/real_dx    
hwidth = real_width/real_dx/2.
length = real_length/real_dx


--- calculate time steps 

rt_time = 2 * length * nref_si

ncyc        = round(rt_time/dt)       -- number of time steps

--- setup source / diagnostic planes for ffts

infty = 10000000000

pulse = { 
   shape   = "Gaussian",
   width   = 100,
   offset  = 0,                    
   attack  = 500,        
   sustain = 0, 
   decay   = 500,   
}

--- excitation pulse

--pulsehwhm   = 75
--pulsehsteps = 500

--- pml cells

cpml = 11      -- in real units it will be real_dx times cpml 

--- setup source / diagnostic planes for ffts

kil  = cpml+10        -- its kiL
kib  = cpml + 5

ki1  = round(length_i1-10)   -- ki1 its k i and number 1 
ki2  = round(length_i-10)
kc  = round(length_i + length_it + length_c/2.)
ko2  = round(length-length_o+10)
ko1  = round(length-cpml-10)

---source injection position definitions

iinj  = round(hwidth_i)             -- x position of source injection
jinj1 = round(height_bsio2)         -- y position 
jinj2 = round(height_bsio2 + height_bsi+height_wg)
kinj  = round(length_i1/10)

ib = round(iinj/2)    -- x position of back reflection

jc           = round(height_bsio2 + (height_bsi+height_wg)/2)
jc1          = round(height_bsio2 )
jc2          = round(height_bsio2 + height_bsi)
jc3          = round(height_bsio2 + height_bsi + height_wg)
js= 2

ici = round (hwidth_i + width_gi + width_ii + width_edw ) 
icc = round (hwidth_c + width_gc + width_cc + width_edw ) 
ico = round (hwidth_o + width_go + width_oo + width_edw ) 

ich = round(hwidth_c/2)   -- center of channel (half half of channel width)

isi = 2
isc = 2
iso = 2

--- contour of waveguide prism (x,y)

wg = {}

wg[1] = {  -30, -30 }
wg[2] = {  hwidth_i , -30} 
wg[3] = {  hwidth_i ,length_i } 
wg[4] = {  hwidth_c ,length_i + length_it }
wg[5] = {  hwidth_c ,length_i + length_it + length_c }
wg[6] = {  hwidth_o ,length_i + length_it + length_c + length_ot }
wg[7] = {  hwidth_o ,length }
wg[8] = {  -30, length }


--- contour of the metal slit


--- print parameters

print("length_i1    = ",  length_i1)
print("length_i2    = ",  length_i2)
print("length_it    = ",  length_it)
print("length_c     = ",  length_c)  
print("length_ot    = ",  length_ot)
print("length_o2    = ",  length_o2)
print("length_o1    = ",  length_o1)
print("hwidth_i     = ",  hwidth_i)   
print("hwidth_c     = ",  hwidth_c)
print("hwidth_o     = ",  hwidth_o)
print("width_g      = ",  width_g) 
print("height_g     = ",  height_g)
print("width_pad    = ",  width_pad)

print("alpha        = ", alpha*180/math.pi)
print("beta         = ", beta*180/math.pi)

print("height_wg    = ",  height_wg)
print("height_bsi   = ",  height_bsi)
print("height_pad   = ",  height_pad)

print("length_i     = ",  length_i)
print("length_o     = ",  length_o)
print("length_t     = ",  length_t)

print("height       = ",  height)    
print("hwidth       = ",  hwidth)
print("length       = ",  length)

print("kil          = ", kib)
print("kib          = ", kil)
print("ki1          = ", ki1)
print("ki2          = ", ki2)
print("kc           = ", kc)
print("ko2          = ", ko2)
print("ko1          = ", ko1)
print("iinj         = ", iinj)
print("jinj1         = ", jinj1)
print("jinj2         = ", jinj2)
print("kinj         = ", kinj)
print("jc           = ", jc)
print("jc1          = ", jc1)
print("jc2          = ", jc2)
print("jc3          = ", jc3)
print("ici          = ", ici)
print("icc          = ", icc)
print("ico          = ", ico)

for i,v in ipairs(wg) do
   print("wg["..tostring(i).."] = ",v[1],v[2])
end

