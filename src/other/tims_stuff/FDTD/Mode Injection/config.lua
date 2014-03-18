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

comp	=	{}
comp.d 	=	{}
comp.t 	=	{}
comp.drude 	=	{}
comp.S 	=	{}

filling 	=	{}
filling.drude 	=	{}

pmls 	=	{}
ebox 	=	{}
Rough 	=	{}

Wire    = 	{}

-- Free space wavelength

real.pump.omega 	=	1.5e15
real.probe.lambda 	= 	1.535 	-- (um)

real.pump.freq 	=	pump.omega/(2*math.pi)
real.pump.lambda 	=	c_SI/pump.freq

real.probe.freq 	=	c_SI/probe.lambda
real.probe.omega 	=	probe.freq*2*math.pi


-- Spatial Resolution

dx 	= 	0.01	-- (um)
w_res 	=	20 	-- Wavelength resolution
res_c 	=	1 	-- Choose between set resolution (1) or wavelength res (2)


-- Simulation Time (s)

real.t.a 	= 	2e-13 		-- Pulse attack time (s)
real.t.s 	= 	2e-11  		-- Pulse sustain time (s)
real.t.d	= 	0	 	-- Pulse decay time (s)
real.t.o	= 	0 		-- Run time after pulse switch off (s)
real.t.w 	=	real.t.a/10 	-- FWHM of pulse (s)


-- Geometry

ntotal 		= 	3 	-- Number of Layers

-- Width of Simulation Box

real.width 	= 	10 	-- (um)

-- Define Layer Thicknesses

real.d.prism  	= 	0.5 	-- (um)
real.d.gap    	= 	0.5 	-- (um)
real.d.layer1 	= 	0.25 	-- (um)
real.d.layer2 	= 	0.47 	-- (um)
real.d.layer3 	= 	1.0 	-- (um)


-- Define Refractive indicies

n 		=	{}
n.bg 		= 	1
n.prism  	= 	1 
n.gap    	= 	1
n.layer1 	= 	1.44
n.layer2 	= 	(11.68)^0.5
n.layer3 	= 	2
n.wire 		=	(9.065)^2

-- Drude Parameters

real.drude.wp   = 	3.13e15 	-- Plasma frequency (rad/s)
real.drude.ga   = 	1.07e14		-- Damping frequency (rad/s)

filling.drude.prism 	= 	0 	-- Filling factors for each layer
filling.drude.gap	=	0
filling.drude.layer1	=	0
filling.drude.layer2 	=	0
filling.drude.layer3	=	0


-- Surface Roughness

Rough.run 	=	0	-- Run surface roughness code
Rough.y 	= 	1 	-- Height variation +/- (nm)
Rough.x 	=	10	-- Width between noise (nm)
Rough.xmin 	=	1	-- Min x variation (nm)
Rough.xmax  	=	2 	-- Max x variation (nm)
Rough.ytype 	=	2 	-- Y displacement type (1 = none, 2 = uniform, 3 = clipped gaussian)
Rough.xtype 	=	1 	-- X displacement type (1 = none, 2 = uniform, 3 = clipped gaussian)
Rough.inter 	=	'all' 	-- Which interfaces to apply noise to

-- Wire Parameters

wire.num 	=	5 	-- Number of wire
wire.h 		=	0.02
wire.w 		=	0.1
wire.s		=	0.1
wire.xst 	=	1


-- Source

amp 		= 	5e8
comp.S.num	=	1	-- Number of sources
real.S.p2p	=	30 	-- Distance between sources (lambda0)
inc_ang 	= 	0 	-- Incident angle (Degrees)

real.S.xst 	= 	1 						-- Source start x position
real.S.length 	= 	0						-- Source length in x direction
real.S.yst 	= 	0.1-real.d.layer1-real.d.layer2-real.d.layer3	-- Source start y position
real.S.yfi 	= 	real.d.gap 		 			-- Source end y position
real.S.lambda 	= 	pump.lambda					-- Permitivity of surrounding material


src_out_f 	=	"source_profile" 		-- Source output filename


-- PML parameters

pmls.t 		= 	11 	-- Thickness of PML Layers
pmls.ext 	= 	2 	-- Extend Layers into PML


----------- Calculated parameters -----------

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


-- Calculate spatial resolution from wavelength res if res_c == 2

n_max 	=	max(n) 	-- Find maximum refractive index

if res_c == 2 then
	dx = real.wavelength/n_max/w_res
end


-- Calculate temporal resolution

dt 		= 	round((0.99 * 1 / math.sqrt(dim)),2) 	-- Courant value
real.dt 	=	dt*dx*1e-6/c_SI				-- Real time step (s)


-- Simulation parameters

sim_t 		=	add(real.t)-real.t.w	-- Total simulation time (s)
ncycles 	=	round(sim_t/real.dt) 	-- Number of cycles
halfcycles 	=	round(ncycles/2) -- Half number of cycles
nstart = 0
fncycles = ncycles
if half.t == 1 then
	if half.v == 1 then
		nstart = 0
		ncycles = halfcycles
		-- Save/load state

		cfg:CHECKPOINT{
				save = true,
				load = false,
				detail = 3,
		}
	else
		nstart = halfcycles+1
		-- Save/load state

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

pulse_off 	=	sim_t-real.t.o
pulse_off 	=	round(pulse_off/real.dt)
pulse_sus 	=	round(real.t.a/real.dt)


-- Convert to computational units

comp.probe.lambda 	=	real.probe.lambda/dx
comp.pump.lambda 	=	real.pump.lambda/dx
comp.width 		=	round(real.width/dx)

-- Geometry

for k,v in pairs(real.d) do
   	comp.d[k] = (v/dx)
end

for k,v in pairs(real.wire) do
   	comp.wire[k] = (v/dx)
end

comp.height 	=	round(real.height/dx) 			-- Height of structure in grid cells
comp.hpg 	=	round(comp.d.prism+comp.d.gap)		-- Height of prism and gap in grid cells

-- Source

-- Calculate sigma from FWHM in terms of wavelength
real.S.sigma 	=	(real.S.sigma*real.probe.lambda)/(2*math.sqrt(2*math.log(2)))

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

real.S.p2p 	=	real.S.p2p*real.probe.lambda
comp.S.p2p 	=	round(real.S.p2p/dx)
C_xst 	= 	{}
C_xst[1] 	= 	comp.S.xst

for n=2, comp.S.num, 1 do
	C_xst[n] = C_xst[n-1]+comp.S.p2p
end
print('width = ',comp.width)
--comp.width 	=	C_xst[comp.S.num]+comp.S.length+C_xst[1]

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

-- Create wire boxes

Bx_wire = {}


for n=1, wire.num, 1 do
	a = string.format("%s%d","wire",n)
	wirestx = comp.wire.xst+(comp.wire.s+comp.wire.w)*(n-1)
	wirefix = wirestx+comp.wire.w
	Bx_wire[a] = Box{
	    		from={wirestx, 0, -infty},
	    		to={wirefix, comp.wire.h, infty}
	    		}
end



-- Create Geometries

scene = Scene{ value = n.bg^2 }

for k,v in pairs(comp.d) do
	
  scene:add{ Bx[k], depth = 1, value = n[k]^2 }

end

for n=1, wire.num, 1 do
        a = string.format("%s%d","wire",n)
	scene:add{Bx_wire[a], depth=0, value = n[wire]^2}
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

a = string.format("wire")
grid[a] = Grid{
	       from={ir0[1], -5,0},
	       to={ir0[2], comp.wire.h+5,0}
	        }

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
a = string.format("Wire")
cfg:CREATE_GEO{
          a,
          scene=scene,
          grid=grid[a],
          }


grid_prev = Grid{
		from={ir0[1], jr0[1],0},
		to={ir0[2], jr0[2],0}
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

for n=1, wire.num, 1, do
        a = string.format("%s%d","wire",n)
	scene_filling:add{Bx_wire[a], depth=0, value = 1}
end

a = string.format("Wire")
cfg:CREATE_GEO{
          a,
          scene=scene_filling,
          grid=grid[a],
          }

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

-- Define energy boxes

ebox.st1 	=	ib[1]+ebox.buf
ebox.fi1 	=	round(comp.width/2)-round(ebox.width/2)+ebox.pos
ebox.st2 	=	ebox.fi1
ebox.fi2 	=	ebox.st2+ebox.width
ebox.st3 	=	ebox.fi2
ebox.fi3 	=	ib[2]-ebox.buf


--- ASSEMBLE CONFIG

cfg:GRID{
   dim = dim,                       -- number of dimensions
   ncyc = { nstart, ncycles },                  -- number of time steps
   dt = dt,                         -- time step in units of cell size
   irange = ir,                     -- i range 
   jrange = jr,                     -- j range  
   krange = kr,
}
print("nstart 	= ",nstart)
print("ncycles 	= ",ncycles)
cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },  -- 0: pec, 1: pml, 2: pmc, 3: pbc

   PML{
      cells = pmls.t,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
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

cfg:FDTD{

   EPSILON{
      REG{
	BOX( box ),
	LOAD_GEO( load )
      }
   },

 OUT{
      file = { "VTK", "Hz" },
      type = { "Hz", "N" },
      time = {1, ncycles, round(ncycles/20)},
      REG{
         BOX{
            { ir[1], ir[2], 2, jr[1], jr[2], 2 }
         }
      },
      on = true
   },
--[[    OUT{
      file = { "GPL", "Hz" },
      type = { "Hz", "N" },
      time = { 21414, ncycles, 342 },
      REG{
         BOX{
            { ir[1], ir[2], 1, mid, mid, 1 }
         }
      },
      on = true
   },
   OUT{
      file = { "GPL", "Point_Hz" },
      type = { "Hz", "N", ".F." },
      time = { 21414, ncycles, 4 },
      REG{
         BOX{
            { midx, midx, 1, mid, mid, 1 }
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


-- Source
for n=1, S_num, 1 do
	name = string.format("source_profile%s.txt",n)
	file_n = io.output(name)
	C_S.xst = C_xst[n]
	for k,v in pairs(S_n) do
 		a = string.format("%s ",v)
		file_n:write(a,C_S[v],"\n")
	end
end

-- Run matlab code
fh = io.open("run.m","w")
cmd = "Field( "..tostring(real.pump.omega)..","
   ..tostring(C_S.yst)..","..tostring(C_S.yfi)..","..tostring(C_S.xst)..","
   ..tostring(dx)..");"
fh:write(cmd.."\n")
cmd = "quit;"
fh:write(cmd.."\n")
fh:close()

os.execute("matlab -r run -nosplash -nodesktop")
fh = io.open("Neff.txt","r")
neff = fh:read("*n")
print("neff = ",neff)
fh = io.open("run2.m","w")
cmd = "Mode( "..tostring(dx)..","
   ..tostring(dt)..");"
fh:write(cmd.."\n")
cmd = "quit;"
fh:write(cmd.."\n")
fh:close()

name = string.format("source_prism_data")
print("name = ",name)
cfg:SRC{
   TFSFINJ{
      invlambda = 1/comp.pump.lambda,
      amplitude = amp*(math.sqrt(8.854187817e-12)*(dx*1e-6)^(3/2)/c_SI),
      pulse = { 
         shape="Gaussian",
         width=comp_pulse_width,
         offset=0,
         attack=comp_attack,
         sustain=0,
         decay=comp_attack
      },
      planewave = { phi=0, theta=90, psi=0, nrefr=neff }
   },
   REG{
      LOAD{ name }
   },
   on = true
}


-- Drude

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

print("resolution 	= ",w_res)
print("Dx 		= ",dx)
print("Dt 		= ",real.dt)
print("Courant 	= ",dt)
print("Lambda 		= ",real.wavelength)
print("Invlambda 	= ",1/comp.wavelength)
print("Height 		= ",jb[2]-jb[1])
print("Width 		= ",ib[2]-ib[1])
print("Height+PMLs 	= ",jr[2]-jr[1])
print("Width+PMLs 	= ",ir[2]-ir[1])
print("Total size 	= ",(jr[2]-jr[1])*(ir[2]-ir[1]))
print("Injection plane 	= ",comp.S.yst)

fh 	= 	io.open("sim_params.txt","w")
if res_c == 2 then
	cmd 	=	"res "..tostring(w_res)
else 
	cmd 	=	"dx "..tostring(dx)
end
fh:write(cmd.."\n")
cmd 	=	"xstart "..tostring(pmls.t)
fh:write(cmd.."\n")
cmd 	=	"xend "..tostring(ir[2]-ir[1]-pmls.t)
fh:write(cmd.."\n")
cmd 	=	"lambda "..tostring(real.wavelength)
fh:write(cmd.."\n")
fh:close()
top 	=	round(comp.inter.prism)+1

-- Reflection box

-- CREATE: config.<part>.in
cfg:CREATE()
