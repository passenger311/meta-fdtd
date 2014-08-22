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

-- Sum Function

function size(a)
	local sum = 0
	for k,v in pairs(a) do
		sum = sum + 1
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
half.t = 1
run = 1

-- Initialise structures

real 	=	{} 		
real.d 	=	{}
real.t 	=	{}
real.drude 	=	{}
real.S 	=	{}
real.pump 	=	{}
real.probe 	=	{}
real.flvl	=	{}

comp	=	{}
comp.d 	=	{}
comp.t 	=	{}
comp.drude 	=	{}
comp.S 	=	{}
comp.pump 	=	{}
comp.probe 	=	{}
comp.flvl 	=	{}

filling 	=	{}
filling.drude 	=	{}

pmls 	=	{}
ebox 	=	{}
Rough 	=	{}


-- Spatial Resolution

dx 	= 	0.01	-- (um)
w_res 	=	20 	-- Wavelength resolution
res_c 	=	1 	-- Choose between set resolution (1) or wavelength res (2)


-- Simulation Time (s)

real.t.s 	= 	10e-12 		-- Pulse sustain time (s)

-- Ignore these for now
real.t.a 	= 	0 		-- Pulse attack time (s)
real.t.d	= 	0	 	-- Pulse decay time (s)
real.t.o	= 	0 		-- Run time after pulse switch off (s)
real.t.w 	=	0 		-- FWHM of pulse (s)


-- Geometry

ntotal 		= 	5 	-- Number of Layers

-- Width of Simulation Box

real.width 	= 	12	-- (um)

-- Define Layer Thicknesses

real.d.layer1	=	0.2
real.d.layer2  	= 	0.1 	-- (um)
real.d.layer3  	= 	0.5 	-- (um)
real.d.layer4 	= 	0.29 	-- (um)
real.d.layer5 	= 	0.6 	-- (um)


-- Define Refractive indicies

n 		=	{}
n.bg 		= 	1
n.layer1        =       1
n.layer2 	=	(11.68)^0.5
n.layer3 	= 	2
n.layer4 	= 	(11.68)^0.5
n.layer5 	= 	2


-- Drude Parameters

real.drude.wp   = 	3.13e15 	-- Plasma frequency (rad/s)
real.drude.ga   = 	1.07e14/10		-- Damping frequency (rad/s)

filling.drude.layer1    =       0       -- Filling factors for each layer
filling.drude.layer2 	= 	0
filling.drude.layer3	=	1
filling.drude.layer4	=	0
filling.drude.layer5	=	1


-- 4lvl Parameters

real.flvl.width = 	1000/1000		-- Width of gain region (um)
real.flvl.layer = 	4 		-- Layer to put gain in
real.flvl.buf 	=	1 		-- Buffer between edge of gain and layer interface (grid cells)

real.flvl.wE 	=	1.55		-- Emission wavelength (um)
real.flvl.wA 	=	1.345 		-- Absorption wavelength (um) Probably dont need to change as it is electrically pumped
real.flvl.wFT 	=	real.flvl.wE 	-- Diag fourier transform wavelength (um)
real.flvl.rp 	=	1e12 		-- Pump rate (1/s)
real.flvl.tuaA = 	40e-15		-- Inverse absorption resonance half width (s)
real.flvl.tuaE = 	40e-15		-- Inverse emission resonance half width (s)
real.flvl.d12 	=	0.4e-9		-- Emission dipole length (m)
real.flvl.d03 	=	0		-- Absorption dipole length (m)
real.flvl.dens 	=	2e18		-- Dipole density (1/cm^3)
real.flvl.stN1 	=	0 		-- Initial population levels
real.flvl.stN2 	=	0		-- Set == 0 so all carriers are initially in the ground state
real.flvl.stN3 	=	0
real.flvl.tau1 = 	100e-15 	-- Energy level lifetimes (s) 1-0 fast
real.flvl.tau2 = 	500e-12		-- 2-1 slower
real.flvl.tau3 = 	100e-15 	-- 3-2 fast
real.flvl.volfac = 	1e8 		-- Related to langevin noise terms, set large
real.flvl.LFE 	=	0 		-- Local field effect (0-off, 1-on) best to leave this off

out_n 	=	50 			-- Number of output files
stab_time = 15e-12

-- PML parameters

pmls.t 		= 	30 	-- Thickness of PML Layers
pmls.ext 	= 	0 	-- Extend Layers into PML

-- Surface Roughness

Rough.run 	=	0	-- Run surface roughness code
Rough.y 	= 	1 	-- Height variation +/- (nm)
Rough.x 	=	10	-- Width between noise (nm)
Rough.xmin 	=	1	-- Min x variation (nm)
Rough.xmax  	=	2 	-- Max x variation (nm)
Rough.ytype 	=	2 	-- Y displacement type (1 = none, 2 = uniform, 3 = clipped gaussian)
Rough.xtype 	=	1 	-- X displacement type (1 = none, 2 = uniform, 3 = clipped gaussian)
Rough.inter 	=	'all' 	-- Which interfaces to apply noise to


----------- Calculated parameters -----------


-- Geometry

ntotal = size(real.d) 			-- Number of layers
real.height 	= 	add(real.d) 	-- Calculate structure height

-- Layer names

layers 		=	{}

for loop=1, ntotal, 1 do
	a	= string.format("%s%d","layer",loop)
	layers[a] = loop+2 
	layers[loop+2] = a
end


-- Calculate spatial resolution from wavelength res if res_c == 2

n_max 	=	max(n) 	-- Find maximum refractive index

if res_c == 2 then
	dx = real.wavelength/n_max/w_res
end

-- Calculate temporal resolution

dt 		= 	round((0.99 * 1 / math.sqrt(dim)),2) 	-- Courant value
real.dt 	=	dt*dx*1e-6/c_SI				-- Real time step (s)

-- Drude parameters

real.drude.lam 		=	{}
real.drude.lam.ga 	=	1e6*c_SI/real.drude.ga 			-- Damping coefficient 
real.drude.lam.wp 	=	1e6*c_SI/(real.drude.wp/(2*math.pi)) 	-- Plasma wavelength

------------ Convert to computational units ---------------------------------

-- Read angle from file

if ang_c == 2 then
	a 	=	string.format("Angle.txt")
	file_n 	=	io.input(a)
	inc_ang = 	io.read("*number")
	file_n:close()
end

-- Simulation parameters

sim_t 		=	add(real.t)-real.t.w	-- Total simulation time (s)
ncycles 	=	round(sim_t/real.dt) 	-- Number of cycles
ncyc_meas	=	5160

pulse_off 	=	sim_t-real.t.o
pulse_off 	=	round(pulse_off/real.dt)
pulse_sus 	=	round(real.t.a/real.dt)
stab_time 	=	round(stab_time/real.dt)

nstart = 2 -- (run - 1)*(ncycles + 1)
nfinish = (ncycles * run) + run - 1
if half.t == 1 then
	if run == 1 then
		-- Save/load state
		cfg:CHECKPOINT{
				save = true,
				load = true,
				detail = 3,
		}
	else
		-- Save/load state
		cfg:CHECKPOINT{
				save = true,
				load = true,
				detail = 3,
		}
	end
end


for k,v in pairs(real.t) do
   	comp.t[k] = round(v/real.dt)
end

-- Geometry

comp.width 		=	round(real.width/dx)

for k,v in pairs(real.d) do
   	comp.d[k] = (v/dx)
end

comp.height 	=	round(real.height/dx) 		-- Height of structure in grid cells
comp.hpg 	=	round(comp.d.layer1)		-- Height of top layer


-- Drude parameters

comp.drude.lam 		=	{}
comp.drude.lam.ga 	= 	real.drude.lam.ga/dx
comp.drude.lam.wp 	=	real.drude.lam.wp/dx
comp.drude.inv 		=	{}
comp.drude.inv.ga 	=	1/comp.drude.lam.ga
comp.drude.inv.wp 	=	1/comp.drude.lam.wp


-- Define computational domain

jstart 	= 	comp.hpg-comp.height 			-- Set zero at interface between layer 1 and 2
jend 	= 	comp.hpg

ib 	= 	{ 0, comp.width }			-- Comp domain
jb 	= 	{ jstart, jend }

ir 	= 	{ ib[1]-pmls.t, ib[2]+pmls.t }  	-- Total computational domain with pmls
jr 	= 	{ jb[1]-pmls.t, jb[2]+pmls.t }  		

ir0 	= 	{ ir[1]-pmls.ext, ir[2]+pmls.ext } 	-- Extend comp domain into PMLs
jr0 	= 	{ jr[1]-pmls.ext, jr[2]+pmls.ext }


-- Define Boxes

-- Calculate interface positions

comp.inter0 		=	jr0[2]
comp.inter1 		=	0
int 	=	0

for loop=2, ntotal, 1 do
   	a 	 	= 	string.format("%s%d","layer",loop)
	b 		=	string.format("%s%d","inter",loop)
   	int 	 	= 	int - comp.d[a]
   	comp[b] 	=  	int

end

comp[b] 	=	jr0[1]  	-- Set final interface to bottom of comp domain	
a2 	 	= 	string.format("%s%d","inter",ntotal-1)
comp.d[a] 	=	comp[a2]-comp[b]

-- Create boxes

Bx = {}

for loop=1, ntotal, 1 do

   	a = string.format("%s%d","layer",loop)
   	fi = string.format("%s%d","inter",loop-1)
	st = string.format("%s%d","inter",loop)

   	Bx[a] 	=    Box{
	    		from={ir0[1]-1, comp[st], -infty},
	    		to={ir0[2]+1, comp[fi], infty}
	    		}
end


-- Create Geometries

scene = Scene{ value = n.bg^2 }

for k,v in pairs(comp.d) do
	
  scene:add{ Bx[k], depth = 1, value = n[k]^2 }

end

-- Gridding

grid = {}

for loop=1, ntotal-1, 1 do
	a = string.format("%s%d","inter",loop)
	grid[a] = Grid{
		      from={ir0[1], round(comp[a]-5),0},
                      to={ir0[2], round(comp[a]+5),0}
                      }
end


-- Geometries and previews

for loop=1, ntotal-1, 1 do
  	a = string.format("%s%d","inter",loop)
   	cfg:CREATE_GEO{
            a,
            scene=scene,
            grid=grid[a],
	    comps=3,
            }

      	cfg:CREATE_PREVIEW{
            a,
            scene=scene,
            grid=grid[a],
            on=true
            }
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


for loop=1, ntotal-1, 1 do
  	a = string.format("%s%d","filling_inter",loop)
	b = string.format("%s%d","inter",loop)
   	cfg:CREATE_GEO{
            a,
            scene=scene_filling,
            grid=grid[b],
	    comps=3,
            }

      	cfg:CREATE_PREVIEW{
            a,
            scene=scene_filling,
            grid=grid[b],
            on=true
            }
end



-- Grid gain

comp.flvl.width = round(real.flvl.width/dx)

c_width = round(comp.width/2)
g_st 	= round(c_width-(comp.flvl.width/2))
g_fi 	= g_st+comp.flvl.width

a 		=	string.format("%s%d","inter",real.flvl.layer-1)
b 		=	string.format("%s%d","inter",real.flvl.layer)
g_yst = comp[b]+real.flvl.buf
g_yfi = comp[a]-real.flvl.buf
g_meas = math.floor(g_yst/5+g_yfi*4/5+0.5)
   Bx.gain = Box{
	    from={g_st, g_yst, -infty},
	    to={g_fi, g_yfi, infty}
	    }

scene_gain = Scene{value=0}
scene_gain:add{Bx.gain,depth=1,value=1}

grid_gain = Grid {
   from={g_st,g_yst,0}, to={g_fi,g_yfi,0}
}
cfg:CREATE_GEO{
   "geo_gain",
   scene=scene_gain,
   grid=grid_gain,
   comps=3,
}
grid_gain_prev = Grid {
   from={ib[1],jb[1],0}, to={ib[2],jb[2],0},
--   yee=false
}
cfg:CREATE_GEO{
   "gain_box",
   scene=scene_gain,
   grid=grid_gain_prev,
   comps=3,
}
os.rename('geo_gain_box.in','gain_box.set')


--- ASSEMBLE CONFIG

cfg:GRID{
   dim = dim,                       -- number of dimensions
   ncyc = { nstart, nfinish }, 	    -- number of time steps
   dt = dt,                         -- time step in units of cell size
   irange = ir,                     -- i range 
   jrange = jr,                     -- j range  
   krange = kr,
   dx = { dx*1e-6*1e-40, 1, 1, 1 }
}

cfg:BOUND{
--   config = { 4, 4, 4, 4, 4, 4 },  -- 0: pec, 1: pml, 2: pmc, 3: pbc
   config = { 1, 1, 1, 1, 1, 1 },

   PML{
      cells = pmls.t,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
--      alpha = 2,
--      alphapot = 3.2,
   }
}

box = {}
load = {}

for loop=1, ntotal, 1 do
   	a = string.format("%s%d","layer",loop)
	fi = string.format("%s%d","inter",loop-1)
	st = string.format("%s%d","inter",loop)
	table.insert( box, {ir0[1], ir0[2], 1, round(comp[st]), round(comp[fi]), 1,
		      ":",n[a]^2,n[a]^2,n[a]^2 })
end

for loop=1, ntotal-1, 1 do
	a = string.format("%s%d","inter",loop)
	table.insert( load, a)
end



midy = -round((comp.d.layer3+(comp.d.layer4/2)))
midx = round((ib[1]+ib[2])/2)
leftx = midx - round(comp.flvl.width/2)
topy = -(comp.d.layer3+2)
thirdy = -round((comp.d.layer3+(comp.d.layer4/3)))

cfg:FDTD{

   EPSILON{
      REG{
	BOX( box ),
	LOAD_GEO( load )
      }
   },
   OUT{
      file = { "GPL", "en_sum_full" },
      type = { "En", "S", ".F."},
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { ib[1], ib[2], 1, jb[1], jb[2], 1, 0, 0, 1 }
         }
      },
      on = true,
   },
   OUT{
      file = { "GPL", "en_sum_strip" },
      type = { "En", "S", ".F."},
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
         }
      },
      on = true,
   },
   OUT{
      file = { "GPL", "enH_sum_strip" },
      type = { "EnH", "S", ".F."},
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
         }
      },
      on = true,
   },
   OUT{
      file = { "GPL", "enE_sum_strip" },
      type = { "EnE", "S", ".F."},
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
         }
      },
      on = true,
   },

   OUT{
      file = { "GPL", "en_sum" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         LOAD_GEO{ "geo_gain" }
      },
      on = true
   },
   OUT{
      file = { "GPL", "enH_sum" },
      type = { "EnH", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         LOAD_GEO{ "geo_gain" }
      },
      on = true
   },   OUT{
      file = { "GPL", "enE_sum" },
      type = { "EnE", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         LOAD_GEO{ "geo_gain" }
      },
      on = true
   },
   OUT{
      file = { "GPL", "en_point" },
      type = { "En", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "enH_point" },
      type = { "EnH", "S", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "enE_point" },
      type = { "EnE", "S", ".F." },
      time = { 0, ncycles, 10 },
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
      time = { 0, ncycles, 10 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
   },
   OUT{
      file = { "GPL", "H_point" },
      type = { "H", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         POINT{
            { c_width, g_meas, 0 }
         }
      },
   },
   OUT{
      file = { "VTK", "E2" },
      type = { "E", "N" },
      time = { 0, ncycles, 25000 },
      REG{
         BOX{
            { ir[1], ir[2], 4, jr[1], jr[2], 4, 0, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "VTK", "E" },
      type = { "E", "N" },
      time = { ncycles-ncyc_meas, ncycles, 20 },
      REG{
         BOX{
            { ir[1], ir[2], 2, jr[1], jr[2], 2, 0, 0, 1 }
         }
      },
      on = false
   },
   OUT{
      file = { "VTK", "Ex" },
      type = { "Ex", "N" },
      time = { ncycles-ncyc_meas, ncycles, 20 },
      REG{
         BOX{
            { ir[1], ir[2], 2, g_meas, g_meas, 1, 0, 0, 1 }
         }
      },
      on = true
   },

}



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
		


	--os.execute("matlab -r seed_gen -nosplash -nodesktop") 	-- Create random noise seed file
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


-- Drude

box_d = {}

for loop=1, ntotal-1, 1 do
	a = string.format("%s%d","inter",loop)
	table.insert( load, a)
end

for loop=1, ntotal, 1 do
	a = string.format("%s%d","layer",loop)
	b = string.format("%s%d","inter",loop)
	c = string.format("%s%d","inter",loop-1)

    	if filling.drude[a] == 1 then

    		table.insert( box_d, {ir0[1], ir0[2], 1, round(comp[b]), round(comp[c]), 1,
	     	 	kbeg0, kend0, 1 } )	
	end
end

load_d = {}
for loop=1, ntotal-1, 1 do
	a = string.format("%s%d","filling_inter",loop)
	table.insert( load_d, a)
end




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

-- 4 lvl system

foutput = io.open("invlambda.in","w+")
--foutput:write(dx/real.flvl.wFT,"\n")
foutput:write(dx/2.51372,"\n")
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
      LOAD_GEO{ "geo_gain" }
   },
   OUT{
      file = { "VTK", "dens" },
      type = { "N", "N", ".T." },
      time = { nfinish-5140, nfinish, 20 },
      REG{
         BOX{
            { g_st-5, g_fi+5, 1,
	      comp.inter4-5, comp.inter3+5, 1 }
         }
      },
   },
   OUT{
      file = { "GPL", "dens_pump" },
      type = { "N", "S", ".F." },
      time = { nstart, nfinish, 20 },
      REG{
         BOX{
            { g_st-5, g_fi+5, 1,
	      comp.inter4-5, comp.inter3+5, 1 }
         }
      },
   },
}


samp_fft=2^round(math.log(real.flvl.wE/dx/7)/math.log(2))
cfg:DIAG{
   PSPEC{
      file = "fft_sum",
      time = { nfinish-samp_fft*2^round(math.log(ncycles/samp_fft/2)/math.log(2)), nfinish, samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=90, theta=90, psi=0 }
   },
   REG{
      BOX{
         { c_width-10, c_width+10, 1,
           comp.inter4, comp.inter3, 1 }
      }
   },
   on = true
}

-- Energy Top
cfg:DIAG{
   EBAL {
     time = {nstart, nfinish, 1},
   },
   REG{
      BOX{
         { g_st-1, g_fi+1, 1,
           comp.inter4 - 55, 0, 1 }
      }
   },
   OUT{
      file = { "GPL", "EnI_gain_vertical" },
      type = { "EnI", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           comp.inter4 - 55, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "En_gain_vertical" },
      type = { "En", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           comp.inter4 - 55, 0, 1 }
      }
      },
      on = true
   },
  OUT{
      file = { "GPL", "Ds_gain_vertical" },
      type = { "DS", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           comp.inter4 - 55, 0, 1 }
      }
      },
      on = true
   },
}

-- Energy Strip
cfg:DIAG{
   EBAL {
     time = {nstart, nfinish, 1},
   },
   REG{
      BOX{
         { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
      }
   },
   OUT{
      file = { "GPL", "EnI_strip" },
      type = { "EnI", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
         BOX{
            { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "En_strip" },
      type = { "En", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
         BOX{
            { ib[1], ib[2], 1, comp.inter4, comp.inter3, 1, 0, 0, 1 }
         }
      },
      on = true
   },
  OUT{
      file = { "GPL", "Ds_strip" },
      type = { "DS", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           comp.inter4 - 55, 0, 1 }
      }
      },
      on = true
   },
}

-- Energy Total
cfg:DIAG{
   EBAL {
     time = {nstart, nfinish, 1},
   },
   REG{
      BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 0, 1 }
      }
   },
   OUT{
      file = { "GPL", "EnI_total" },
      type = { "EnI", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 0, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "En_total" },
      type = { "En", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 0, 1 }
      }
      },
      on = true
   },
  OUT{
      file = { "GPL", "Ds_total" },
      type = { "DS", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 0, 1 }
      }
      },
      on = true
   },
}

-- Energy Total2
cfg:DIAG{
   EBAL {
     time = {nstart, nfinish, 1},
   },
   REG{
      BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 20, 1 }
      }
   },
   OUT{
      file = { "GPL", "EnI_total2" },
      type = { "EnI", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 20, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "En_total2" },
      type = { "En", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 20, 1 }
      }
      },
      on = true
   },
  OUT{
      file = { "GPL", "Ds_total2" },
      type = { "DS", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-2000, g_fi+2000, 1,
           comp.inter4 - 55, 20, 1 }
      }
      },
      on = true
   },
}
-- Just Gain
cfg:DIAG{
   EBAL {
     time = {nstart, nfinish, 1},
   },
   REG{
      LOAD_GEO{ "geo_gain" }
   },
   OUT{
      file = { "GPL", "EnI_just_gain" },
      type = { "EnI", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           g_yst-1, g_yfi+1, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "En_just_gain" },
      type = { "En", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           g_yst-1, g_yfi+1, 1 }
         }
      },
      on = true
   },
  OUT{
      file = { "GPL", "Ds_just_gain" },
      type = { "DS", "N", ".F." },
      time = {nstart, nfinish, 10},
      REG{
      	BOX{
         { g_st-1, g_fi+1, 1,
           g_yst-1, g_yfi+1, 1 }
         }
      },
      on = true
   },
}


-- CREATE: config.<part>.in
cfg:CREATE()
