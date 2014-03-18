cfg = CONFIG{scenes=true}

----------- Functions -----------

-- Rounding Function

function round(num, idp)
  	local mult = 10^(idp or 0)
  	return math.floor(num * mult + 0.5) / mult
end

-- Sum Function

function add(a)
	local sum = 0
	for k,v in pairs(a) do
		sum = sum + v
	end
	return sum
end

-- Find Max Function

function max(a)
	local mx = 0
	for k,v in pairs(a) do
   		if v > mx then
      			mx = v
   		end
	end
	return mx
end

-- Find Min Function

function min(a)
	local mn = 1e8
	for k,v in pairs(a) do
   		if v < mn then
      			mn = v
   		end
	end
	return mn
end



----------- Define Simulation constants -----------


dim 	= 	2		-- Number of Dimensions
c_SI 	= 	299792458	-- Speed of Light (m/s)
infty 	= 	1E8 		-- Define Infinity as a Large Number



----------- 	      Inputs 		----------- 

-- First/second half of simulation
half = {}
half.t = 0
half.v = 1

-- Initialise structures

real 	=	{} 		
real.d 	=	{}
real.t 	=	{}
real.drude 	=	{}
real.S 	=	{}
real.flvl	=	{}

comp	=	{}
comp.d 	=	{}
comp.t 	=	{}
comp.drude 	=	{}
comp.S 	=	{}
comp.flvl 	=	{}

filling 	=	{}
filling.drude 	=	{}

pmls 	=	{}
ebox 	=	{}
Rough 	=	{}

-- Number of output files
outNum 	=	200

-- Free space wavelength

real.wavelength 	= 	1.55 	-- (um)


-- Spatial Resolution

dx 	= 	0.01	-- (um)
w_res 	=	20 	-- Wavelength resolution
res_c 	=	1 	-- Choose between set resolution (1) or wavelength res (2)


-- Geometry

ntotal 		= 	3 	-- Number of Layers


-- Width of Simulation Box

real.width 	= 	200 	-- (um)
WinH 		=	20 	-- Width in terms of stack height
w_c		=	1 	-- Choose between set width (1) or height related (2)


-- Define Layer Thicknesses

real.d.prism  	= 	0.4 	-- (um)
real.d.gap    	= 	0.1 	-- (um)
real.d.layer1 	= 	0.5 	-- (um)
real.d.layer2 	= 	0.29 	-- (um)
real.d.layer3 	= 	0.7 	-- (um)


-- Define Refractive indicies

n 		=	{}
n.bg 		= 	1
n.prism  	= 	1 
n.gap    	= 	1
n.layer1 	= 	2
n.layer2 	= 	(11.68)^0.5
n.layer3 	= 	2


-- Drude Parameters

real.drude.wp   = 	3.13e15 	-- Plasma frequency (rad/s)
real.drude.ga   = 	1.07e14		-- Damping frequency (rad/s)

filling.drude.prism 	= 	0 	-- Filling factors for each layer
filling.drude.gap	=	0
filling.drude.layer1	=	1
filling.drude.layer2 	=	0
filling.drude.layer3	=	1

-- Gain Dimensions 
real.flvl.width = 	0.4		-- Width of gain region (um)
real.flvl.layer = 	4 		-- Layer to put gain in
real.flvl.buf 	=	1 		-- Buffer between edge of gain and layer interface (grid cells)


-- 4lvl Parameters
real.flvl.wE 	=	1.55		-- Emission wavelength (um)
real.flvl.wA 	=	1.345 		-- Absorption wavelength (um) Note: unused as electrically pumped
real.flvl.wFT 	=	1.5467 		-- Diag fourier transform wavelength (um)
real.flvl.rp 	=	2e12 		-- Pump rate (1/s)
real.flvl.tuaA = 	40e-15		-- Inverse absorption resonance half width (s)
real.flvl.tuaE = 	40e-15		-- Inverse emission resonance half width (s)
real.flvl.d12 	=	0.4e-9		-- Emission dipole length (m)
real.flvl.d03 	=	0		-- Absorption dipole length (m)
real.flvl.dens 	=	2e18		-- Dipole density (1/cm^3)
real.flvl.stN1 	=	0 		-- Initial population levels
real.flvl.stN2 	=	0.7		-- Set == 0 so all carriers are initially in the ground state
real.flvl.stN3 	=	0
real.flvl.tau1 = 	100e-15 	-- Energy level lifetimes (s) 1-0 fast
real.flvl.tau2 = 	500e-12		-- 2-1 slower
real.flvl.tau3 = 	100e-15 	-- 3-2 fast
real.flvl.volfac = 	1e8 		-- Related to langevin noise terms, set large
real.flvl.LFE 	=	0 		-- Local field effect (0-off, 1-on) best to leave this off


-- Surface Roughness

Rough.run 	=	1	-- Run surface roughness code
Rough.y 	= 	1 	-- Height variation +/- (nm)
Rough.x 	=	10	-- Width between noise (nm)
Rough.xmin 	=	1	-- Min x variation (nm)
Rough.xmax  	=	2 	-- Max x variation (nm)
Rough.ytype 	=	2 	-- Y displacement type (1 = none, 2 = uniform, 3 = clipped gaussian, 4 = sin wave)
Rough.xtype 	=	1 	-- X displacement type (1 = none, 2 = uniform, 3 = clipped gaussian)
Rough.inter 	=	'all' 	-- Which interfaces to apply noise to
--[[roughWavelength =  	5e-6    -- Wavelength of sin wave
roughWavenumber =  	2*math.pi/roughWavelength

a 	= 	string.format("wavenumber.txt")
file_n  = 	io.output(a)
file_n:write(roughWavenumber,"\n")
file_n:close()]]--

a       =       string.format("wavenumber.txt")
file_n  =       io.input(a)
roughWavenumber =       io.read("*number")
file_n:close()



-- Source Parameters Time (s)

real.amp 	= 	1e6
real.t.a 	= 	2e-13 				-- Pulse attack time (s)
real.t.s 	= 	6e-13  				-- Pulse sustain time (s)
real.t.d	= 	real.t.a 			-- Pulse decay time (s)
real.t.o	= 	3e-13 				-- Run time after pulse switch off (s)
real.t.w 	=	real.t.a/10 			-- FWHM of pulse (s)

comp.S.num	=	1				-- Number of sources
real.S.p2p	=	30 				-- Distance between sources (lambda0)
ang_c 		=	2				-- Set angle below (1) or read from file (2)
inc_ang 	= 	72.486 				-- Incident angle (Degrees)


real.S.xst 	= 	1				-- Source start x position
real.S.length 	= 	real.width-(2*real.S.xst)	-- Source length in x direction
real.S.yst 	= 	real.d.gap+0.1			-- Source start y position
real.S.yfi 	= 	real.S.yst 		 	-- Source end y position
real.S.eps 	= 	n.prism^2			-- Permitivity of surrounding material
real.S.mu 	= 	1 				-- Permeability of surrounding material
real.S.sigma 	=	30 				-- FWHM in mulitples of free space wavelength

src_out_f 	=	"source_profile" 		-- Source output filename


-- PML parameters

pmls.t 		= 	30 	-- Thickness of PML Layers
pmls.ext 	= 	0 	-- Extend Layers into PML


-- Gridding buffer
gridBuf 	=	5


----------- Calculated parameters -----------

-- Read angle from file

if ang_c == 2 then
	a 	=	string.format("Angle.txt")
	file_n 	=	io.input(a)
	inc_ang = 	io.read("*number")
	file_n:close()
end


a 	=	string.format("Mag.txt")
file_n 	=	io.input(a)
Rough.y = 	io.read("*number")
file_n:close()

-- Geometry

real.height 	= 	add(real.d) 	-- Calculate structure height

-- Layer names

layers 		=	{}
layers.prism 	= 	1
layers[1] 	=	"prism"
layers.gap 	=	2
layers[2] 	=	"gap"


for n=1, ntotal, 1 do
	a	= string.format("%s%d","layer",n)
	layers[a] = n+2 
	layers[n+2] = a
end


-- Calculate waveguide width from height if w_c = 2

if w_c == 2 then
	real.width = WinH*real.height
end


-- Calculate spatial resolution from wavelength res if res_c == 2

n_max 	=	max(n) 	-- Find maximum refractive index

if res_c == 2 then
	dx = real.wavelength/n_max/w_res
end


-- Calculate temporal resolution

dt 		= 	round((0.99 * 1 / math.sqrt(dim)),2) 	-- Courant value
real.dt 	=	dt*dx*1e-6/c_SI				-- Real time step (s)


-- Simulation parameters

sim_t 		=	add(real.t) - real.t.w	-- Total simulation time (s)
ncycles 	=	round(sim_t/real.dt) 	-- Number of cycles
halfcycles 	=	round(ncycles/2) -- Half number of cycles
T 		=	1/(c_SI/(real.flvl.wE*1e-6))
compT		=	T/real.dt
meas_inc 	=	math.floor(compT/11)
ncyc_meas	=	meas_inc*258
nstart = 0
fncycles = ncycles
if half.t == 1 then
	if half.v == 1 then

		cfg:CHECKPOINT{
				save = true,
				load = false,
				detail = 3,
		}
	else

		cfg:CHECKPOINT{
				save = false,
				load = true,
				detail = 3,
		}
	end
else
	--[[cfg:CHECKPOINT{
			save = true,
			load = false,
			detail = 3,
	}]]--

end


real.drude.lam 		=	{}
real.drude.lam.ga 	=	1e6*c_SI/real.drude.ga 			-- Damping coefficient wavelength
real.drude.lam.wp 	=	1e6*c_SI/(real.drude.wp/(2*math.pi)) 	-- Plasma wavelength




-- Convert to computational units

comp.wavelength 	=	real.wavelength/dx
comp.width 		=	round(real.width/dx)

-- Geometry

for k,v in pairs(real.d) do
   	comp.d[k] = (v/dx)
end

comp.height 	=	round(real.height/dx) 			-- Height of structure in grid cells
comp.hpg 	=	round(comp.d.prism+comp.d.gap)		-- Height of prism and gap in grid cells

-- Source

-- Calculate sigma from FWHM in terms of wavelength
real.S.sigma 	=	(real.S.sigma*real.wavelength)/(2*math.sqrt(2*math.log(2)))

for k,v in pairs(real.S) do
   	comp.S[k] = round(v/dx)
end
comp.S.yst 	=	comp.S.yst-2
comp.S.yfi 	=	comp.S.yst
comp.S.eps 	= 	real.S.eps
comp.S.mu 	= 	real.S.mu


for k,v in pairs(real.t) do
   	comp.t[k] = round(v/real.dt)
end

-- Calculate start positions of sources

real.S.p2p 	=	real.S.p2p*real.wavelength
comp.S.p2p 	=	round(real.S.p2p/dx)
C_xst 	= 	{}
C_xst[1] 	= 	comp.S.xst

for n=2, comp.S.num, 1 do
	C_xst[n] = C_xst[n-1]+comp.S.p2p
end
print('width = ',comp.width)


-- Drude parameters

comp.drude.lam 		=	{}
comp.drude.lam.ga 	= 	real.drude.lam.ga/dx
comp.drude.lam.wp 	=	real.drude.lam.wp/dx
comp.drude.inv 		=	{}
comp.drude.inv.ga 	=	1/comp.drude.lam.ga
comp.drude.inv.wp 	=	1/comp.drude.lam.wp


-- Define computational domain

jstart 	= 	comp.hpg-comp.height
jend 	= 	comp.hpg

ib 	= 	{ 0, comp.width }			-- Comp domain
jb 	= 	{ jstart, jend }

ir 	= 	{ ib[1]-pmls.t, ib[2]+pmls.t }  	-- Total computational domain with pmls
jr 	= 	{ jb[1]-pmls.t, jb[2]+pmls.t }  		

ir0 	= 	{ ir[1]-pmls.ext, ir[2]+pmls.ext } 	-- Extend comp domain into PMLs
jr0 	= 	{ jr[1]-pmls.ext, jr[2]+pmls.ext }


-- Define Boxes

-- Calculate interface positions

comp.inter 		=	{}
comp.inter.top 		=	jr0[2]
comp.inter.prism 	=	comp.d.gap
comp.inter.gap 		=	comp.inter.prism-comp.d.gap
int 			=	comp.inter.gap
comp.d.prism 		=	comp.inter.top-comp.inter.prism

for n=1, ntotal, 1 do
   	a 	 	= 	string.format("%s%d","layer",n)
   	int 	 	= 	int - comp.d[a]
   	comp.inter[a] 	=  	int

end

comp.inter[a] 	=	jr0[1]  	-- Set final interface to bottom of comp domain	
a2 	 	= 	string.format("%s%d","layer",ntotal-1)
comp.d[a] 	=	comp.inter[a2]-comp.inter[a]

-- Create boxes

Bx = {}

Bx.prism 	=    Box{
             		from={ir0[1], comp.inter.prism, -infty},
             		to={ir0[2], comp.inter.top, infty}
	     	        }

Bx.gap 		=    Box{
      			from={ir0[1], comp.inter.gap, -infty},
  			to={ir0[2], comp.inter.prism, infty}
			}

for n=1, ntotal, 1 do

   	a = string.format("%s%d","layer",n)
   	
   	if n == 1 then
     		b = "gap"
	else
		b = string.format("%s%d","layer",n-1)	
   	end

   	Bx[a] 	=    Box{
	    		from={ir0[1], comp.inter[a], -infty},
	    		to={ir0[2], comp.inter[b], infty}
	    		}
end


-- Create Geometries

scene = Scene{ value = n.bg^2 }

for k,v in pairs(comp.d) do
	
  scene:add{ Bx[k], depth = 1, value = n[k]^2 }

end

-- Gridding

grid = {}

for k,v in pairs(comp.d) do
	a = string.format("%s",k)
	if a ~= string.format("%s%d","layer",ntotal) then
    		grid[k]  =  Grid{
		  		from={ir0[1], round(comp.inter[k])-5,0},
		  		to={ir0[2], round(comp.inter[k])+5,0}
		  		}
	end
end

-- Geometries and previews

for k,v in pairs(comp.d) do
   a = string.format("%s",k)
   if a ~= string.format("%s%d","layer",ntotal) then
      cfg:CREATE_GEO{
          a,
          scene=scene,
          grid=grid[k],
          }

      cfg:CREATE_PREVIEW{
          a,
          scene=scene,
          grid=grid[k],
          on=true
          }
   end
end

grid_prev = Grid{
		from={20, jr0[1],0},
		to={30, jr0[2],0}
		}
cfg:CREATE_PREVIEW{
		   "slice",
	            scene=scene,
		    grid=grid_prev,
		   }

-- Grid filling factors

-- Define filling factors

scene_filling = Scene{depth=1, value=0}

for k,v in pairs(filling.drude) do
	if v == 1 then
   		scene_filling:add{Bx[k],depth = 1, value = 1}
	end
end


nms = {}
ind2 = 0

for k,v in pairs(filling.drude) do

	if v == 1 then
		
		ind 	=	layers[k]
		if ind == ntotal+2 then
			nm 	=	layers[ind-1]
			ind2 	=	ind-1
		elseif ind == 1 then
			nm 	=	layers[ind+1]
			ind2 	=	ind-1
		else
			nm 	=	layers[ind-1]
			ind2 	=	ind-1
			a = string.format("filling_%s_interface",nm)
			nms[ind2] = a
   			cfg:CREATE_GEO{
		    	    a,
		     	    scene=scene_filling,
	            	    grid=grid[nm],
		    	    }
			nm 	=	layers[ind]
			ind2 	=	ind
		end
		
		a = string.format("filling_%s_interface",nm)
		nms[ind2] = a
   		cfg:CREATE_GEO{
		    a,
		    scene=scene_filling,
	            grid=grid[nm],
		    }

	
	end
end 

-- Grid gain

comp.flvl.width = (real.flvl.width/dx)
c_width = 	(comp.width/2)
g_st 	= 	(c_width-(comp.flvl.width/2))
g_fi 	= 	g_st+comp.flvl.width

a	=	string.format("%s%d","layer",real.flvl.layer-1)
b	=	string.format("%s%d","layer",real.flvl.layer)
g_yst 	= 	comp.inter.layer2+real.flvl.buf
g_yfi 	= 	comp.inter.layer1-real.flvl.buf
g_meas 	= 	math.floor(g_yst/5+g_yfi*4/5+0.5)

Bx.gain = Box{
	     from={g_st, g_yst, -infty},
	     to={g_fi, g_yfi, infty}
	     }

scene_gain = Scene{value=0}
scene_gain:add{Bx.gain,depth=1,value=1}

g_st 	=	round(g_st)
g_fi 	=	round(g_fi)
g_yst 	=	round(g_yst)
g_yfi 	=	round(g_yfi)

grid_gain = Grid {
   		 from={g_st-gridBuf,g_yst-gridBuf,0},
		 to={g_fi+gridBuf,g_yfi+gridBuf,0}
		 }

cfg:CREATE_GEO{
   	      "Gain",
   	      scene=scene_gain,
   	      grid=grid_gain,
   	      comps=3,
}

cfg:CREATE_PREVIEW{
            "Gain",
            scene=scene_gain,
            grid=grid_gain,
            on=true
            }




--- ASSEMBLE CONFIG

cfg:GRID{
   dim = dim,                       -- number of dimensions
   ncyc = { nstart, ncycles },                  -- number of time steps
   dt = dt,                         -- time step in units of cell size
   irange = ir,                     -- i range 
   jrange = jr,                     -- j range  
   krange = kr,
   dx = { dx*1e-6*1e-40, 1, 1, 1 }
}
print("nstart 	= ",nstart)
print("ncycles 	= ",ncycles)
--[[cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },  -- 0: pec, 1: pml, 2: pmc, 3: pbc

   PML{
      cells = pmls.t,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}]]--

cfg:BOUND{
   config = { 4, 4, 4, 4, 4, 4 },  -- 0: pec, 1: pml, 2: pmc, 3: pbc
   CPML{
      cells = pmls.t,
      pot = 3.2,
      sigma = 4.94444444444444,
      kappa = 1.1,
      alpha = 2,
      alphapot = 3.2,
   }
}

box = {}
load = {}

for k,v in pairs(comp.d) do
    table.insert( box, {ir0[1], ir0[2], 1, round(comp.inter[k]), round(comp.inter[k]+comp.d[k]), 1,
	     	 ":",n[k]^2,n[k]^2,n[k]^2 } )
    if k ~= string.format("%s%d","layer",ntotal) then
	table.insert( load, string.format("%s",k) )
    end    	
end

mid = -round((comp.d.layer1+(comp.d.layer2/2)))
midx = round((ir[1]+ir[2])/2)

tst = 1.05e-12
tend = 1.1e-12
comptst = round(tst/real.dt)
comptend = round(tend/real.dt)
print('real dt = ',real.dt)
period = real.wavelength*1e-6/299792458
compPeriod = period/real.dt
num = 20
tinc = round(compPeriod/num)
tfi = round(comptst + 2*compPeriod)

cfg:FDTD{

   EPSILON{
      REG{
	BOX( box ),
	LOAD_GEO( load )
      }
   },

--[[  OUT{
      file = { "GPL", "en_point" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "E_point" },
      type = { "E", "N", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
   },
   OUT{
      file = { "GPL", "H_point" },
      type = { "H", "N", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
   },
   OUT{
      file = { "VTK", "En" },
      type = { "En", "N" },
      time = { ncycles-ncyc_meas, ncycles, meas_inc },
      REG{
         BOX{
            { ir[1], ir[2], 2, comp.inter.layer2-20, comp.inter.layer1+20, 1, 0, 0, 1 }
         }
      },
      on = true
   },]]--
   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { comptend, comptend, 1 },
      REG{
         BOX{
            { ir[1], ir[2], 4, jr[1], jr[2], 2, 0, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "VTK", "H" },
      type = { "H", "N" },
      time = { comptend, comptend, 1 },
      REG{
         BOX{
            { ir[1], ir[2], 4, jr[1], jr[2], 2, 0, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "VTK", "S" },
      type = { "S", "N" },
      time = { comptend, comptend, 1 },
      REG{
         BOX{
            { ir[1], ir[2], 4, jr[1], jr[2], 2, 0, 0, 1 }
         }
      },
      on = true
   },
--[[   OUT{
      file = { "VTK", "Ex" },
      type = { "Ex", "N" },
      time = { 0, ncycles, round(ncycles/outNum) },
      REG{
         BOX{
            { ir[1], ir[2], 2, jr[1], jr[2], 2, 0, 0, 1 }
         }
      },
      on = true
   },]]--

}

-- Source

for loop=1, comp.S.num, 1 do

	a 	= 	string.format("%s%s.txt",src_out_f,loop)

	file_n  = 	io.output(a)

	S_n 	= 	{}
	S_n[1] 	= 	"xst"
	S_n[2] 	= 	"length"
	S_n[3] 	= 	"yst"
	S_n[4] 	= 	"yfi"
	S_n[5] 	= 	"eps"
	S_n[6] 	= 	"mu"
	S_n[7] 	= 	"sigma"
	S_n[8] 	=	"num"
	comp.S.xst 	=	C_xst[loop]

	for k,v in pairs(S_n) do
	 	a = string.format("%s ",v)
		file_n:write(a,comp.S[v],"\n")
	end
	
	file_n:close()

end

-- Run Matlab codes

os.execute("matlab -r prismsource -nosplash -nodesktop") 	-- Create prism source file

-- Create surface roughness
if Rough.run ~= 0 then
	a 	=	'R_interfaces.txt'
	file_n 	=	io.output(a)
	if Rough.inter == 'all' then

		for k,v in pairs(comp.d) do
			a = string.format("%s",k)
			if a ~= string.format("%s%d","layer",ntotal) then
				epsa 	=	(n[layers[layers[a]+1]])^2
				epsb 	=	(n[a])^2
				
				file_n:write(a.."\t"..tostring(round(comp.inter[k],10)).."\t"..tostring(layers[k]).."\t"..tostring(epsa).."\t"..tostring(epsb).."\n")	
			end
		end

		for k,v in pairs(nms) do
			a = string.format("%s",v)
			k2 = layers[k]
			epsa 	=	filling.drude[layers[layers[k2]+1]]
			epsb 	=	filling.drude[k2]
			
			file_n:write(a.."\t"..tostring(round(comp.inter[k2],10)).."\t"..tostring(layers[k2]).."\t"..tostring(epsa).."\t"..tostring(epsb).."\n")

		end
	else
		for k,v in pairs(Rough) do
			if k ~= run then
				a = string.format("%s",k)
				if a ~= string.format("%s%d","layer",ntotal) then
					epsa 	=	(n[layers[layers[a]+1]])^2
					epsb 	=	(n[a])^2
					file_n:write(a.."\t"..tostring(comp.inter[k]).."\t"..tostring(layers[k]).."\t"..tostring(epsa).."\t"..tostring(epsb).."\n")
				end
			end
		end

		for k,v in pairs(Rough) do
			if k ~= run then
				k2 = string.format("filling_%s_interface",k)
				epsa 	=	filling.drude[layers[layers[k2]+1]]
				epsb 	=	filling.drude[k2]
				file_n:write(a.."\t"..tostring(comp.inter[k2]).."\t"..tostring(layers[k2]).."\t"..tostring(epsa).."\t"..tostring(epsb).."\n")
			end

		end

	end
	file_n:close()
		


	os.execute("matlab -r seed_gen -nosplash -nodesktop") 	-- Create random noise seed file
	a 	=	'R_params.txt'
	fh 	= 	io.output(a)
	cmd 	=	"Width = "..tostring(ir0[2]-ir0[1])
	fh:write(cmd.."\n")
	cmd 	=	"x noise min = "..tostring(Rough.xmin)
	fh:write(cmd.."\n")
	cmd 	=	"x noise max = "..tostring(Rough.xmax)
	fh:write(cmd.."\n")
	cmd 	=	"y noise = "..tostring(Rough.y)
	fh:write(cmd.."\n")
	cmd 	=	"x type = "..tostring(Rough.xtype)
	fh:write(cmd.."\n")
	cmd 	=	"y type = "..tostring(Rough.ytype)
	fh:write(cmd.."\n")
	cmd 	=	"n_max = "..tostring(n_max)
	fh:write(cmd.."\n")
	cmd 	=	"Lambda = "..tostring(real.wavelength)
	fh:write(cmd.."\n")
	if res_c == 2 then
		cmd 	=	"res = "..tostring(w_res)
	else 
		cmd 	=	"dx = "..tostring(dx)
	end
	fh:write(cmd.."\n")
	cmd 	=	"x width = "..tostring(Rough.x)
	fh:write(cmd.."\n")
	fh:close()
	os.execute("matlab -r Roughness2 -nosplash -nodesktop")	-- Run roughness code
end


box_d = {}

for k,v in pairs(filling.drude) do
    	if v == 1 then
    		table.insert( box_d, {ir0[1], ir0[2], 1, round(comp.inter[k]), round(comp.inter[k]+comp.d[k]), 1,
	     	 	kbeg0, kend0, 1 } )	
	end
end

load_d = {}
for k,v in pairs(nms) do
	a = string.format("%s",nms[k])
	table.insert( load_d, nms[k] )
end


-- Source
comp.wavelength 	=	real.wavelength/dx
comp.amp        =       real.amp*(math.sqrt(8.854187817e-12)*(dx*1e-6)^(3/2)/c_SI)*1e20

for loop=1, comp.S.num, 1 do

	a 	= 	string.format("source_prism_data%s",loop)
	print(a)
	cfg:SRC{
   		TFSFINJ{
      			invlambda = 1/comp.wavelength,
			amplitude = comp.amp,
      			pulse = { 
         			shape="Gaussian",
         			width=comp.t.w,
         			offset=0,
         			attack=comp.t.a,
         			sustain=comp.t.s,
         			decay=comp.t.d
      				},
      			planewave = { phi=-inc_ang, theta=90, psi=0, nrefr=n.prism }
   			},
   		REG{
      		    LOAD{ a }
   		    },
   		on = true
		}
end

-- 4 lvl system

--[[foutput = io.open("invlambda.in","w+")
foutput:write(dx/real.flvl.wFT,"\n")
foutput:close()

cfg:MAT{
   LVFOURLVLEP{
      rp = real.flvl.rp*dx*1e-6/c_SI,
      invlambdal = {dx/real.flvl.wE,dx/real.flvl.wA},
      gammal = {(1/real.flvl.tuaE)*1e-6*dx/c_SI,(1/real.flvl.tuaA)*1e-6*dx/c_SI},
      dipole12 = real.flvl.d12*1e6/dx,
      dipole03 = real.flvl.d12*1e6/dx,
      dens = real.flvl.dens*1e6*(dx*1e-6)^3,
      seed = 421761282,
      start = {real.flvl.stN1,real.flvl.stN2,real.flvl.stN3},
      gamma = {0,(1/real.flvl.tau1)*1e-6*dx/c_SI,	
      (1/real.flvl.tau2)*1e-6*dx/c_SI,
      (1/real.flvl.tau3)*1e-6*dx/c_SI},
      volfac = real.flvl.volfac,
      LFE = real.flvl.LFE,
   },
   REG{
      LOAD_GEO{ "Gain" }
   },
  OUT{
      file = { "VTK", "dens" },
      type = { "N", "N", ".T." },
      time = {0, ncycles, round(ncycles/outNum)},
      REG{
         BOX{
            { g_st-5, g_fi+5, 1,
	      comp.inter.layer2-5, comp.inter.layer1+5, 1 }
         }
      },
   },
   OUT{
      file = { "GPL", "dens_pump" },
      type = { "N", "S", ".F." },
      time = { 0, ncycles, 20 },
      REG{
         BOX{
            { g_st-5, g_fi+5, 1,
	      comp.inter.layer2-5, comp.inter.layer1+5, 1 }
         }
      },
   },
}]]--


cfg:MAT{
   on = true,
     DRUDE{
     invlambdapl 	= 	comp.drude.inv.wp,
     gammapl 		= 	comp.drude.inv.ga,
     order = 1.0,
     },
     REG{
        BOX( box_d ),
	LOAD_GEO( load_d ),
         }
}


--[[samp_fft=2^round(math.log(real.flvl.wE/dx/7)/math.log(2))
--samp_fft=2^14
cfg:DIAG{
   PSPEC{
      file = "fft_sum",
      time = { ncycles-samp_fft*2^15, ncycles, samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=90, theta=90, psi=0 }
   },
   REG{
      BOX{
         { c_width-10, c_width+10, 1,
           comp.inter.layer2, comp.inter.layer1, 1 }
      }
   },
   on = true
}

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "ModePlot",
      type = "F",
      time = { ncycles-samp_fft*2^15, ncycles, samp_fft },
   },
   REG{
      BOX{
         { ib[1], ib[2], 2, jb[1], jb[2], 2, 0, 0, 1, }
      }
   },
   on = true
}]]--

yst 	=	jb[1]+10
yfi 	=	1
xst 	=	ib[1]+10
xfi 	=	ib[2]-10

cfg:DIAG{
	   EBAL {
	     time = {0, ncycles, 1},
	   },
	   REG{
	      BOX{
		 {xst, xfi, 1, yst,yfi,1}
	      }
	   },
           OUT{
	      file = { "GPL", "Centre_of_energy_position" },
	      type = { "UR", "N", ".F." },
	      time = {0, ncycles, 5},
	      REG{
	         BOX{
	  	    {xst, xfi, 1, yst,yfi,1}
	         }
	      },
	      on = true
	   },
	   OUT{
	      file = { "GPL", "Second_order_moment" },
	      type = { "URR", "N", ".F." },
	      time = {0, ncycles, 5},
	      REG{
	         BOX{
	  	    {xst, xfi, 1, yst, yfi,1}
	         }
	      },
	      on = true
	   },
	   OUT{
	      file = { "GPL", "Time_Integrated_energy_density" },
	      type = { "EnI", "N", ".F." },
	      time = {0, ncycles, 5},
	      REG{
	         BOX{
	  	    {xst, xfi, 1, yst, yfi,1}
	         }
	      },
	      on = true
	   },
	   OUT{
	      file = { "GPL", "Energy_density" },
	      type = { "En", "N", ".F." },
	      time = {0, ncycles, 5},
	      REG{
	         BOX{
	  	    {xst, xfi, 1, yst, yfi,1}
	         }
	      },
	      on = true
	   },
	   OUT{
	      file = { "GPL", "Poynting_Vector" },
	      type = { "DS", "N", ".F." },
	      time = {0, ncycles, 5},
	      REG{
	         BOX{
	  	    {xst, xfi, 1, yst, yfi,1}
	         }
	      },
	      on = true
	   }
	}



-- CREATE: config.<part>.in
cfg:CREATE()
