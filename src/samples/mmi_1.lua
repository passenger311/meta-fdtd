
cfg = CONFIG{scenes=true}


--- real device parameters

real_wavelength = 1.550 -- um
real_width_mmi = 2.43 -- um
real_length_mmi = 6.15 -- um
real_width_wg = 0.400 -- um
real_height_bsio2 = 0.7
real_height_wg = 0.220
real_length_wg1 = 1.4
real_hwidth_sep = 0.65


nref_si = 3.47
nref_sio2 = 1.444

eps_si = nref_si^2
eps_sio2 = nref_sio2^2

--- device in natural HL, dx=dy=dz=1.

width = 281   -- i
height = 101  -- j
length = 800  -- k

maxi = math.floor( (width - 1) / 2  )
maxj = math.floor( height - 1 )
maxk = math.floor( length - 1 )

resolution = 40

real_dx = real_wavelength / nref_si / resolution

grid_height_wg = math.floor(real_height_wg/real_dx+0.5)
grid_height_bsio2 = math.floor(real_height_bsio2/real_dx)
grid_hwidth_wg = math.floor(real_width_wg/2/real_dx+0.5)
grid_hwidth_mmi = math.floor(real_width_mmi/2/real_dx+0.5)
grid_length_wg1 = math.floor(real_length_wg1/real_dx+0.5)
grid_length_mmi = math.floor(real_length_mmi/real_dx+0.5)
grid_hwidth_sep = math.floor(real_hwidth_sep/real_dx+0.5)

print("irange (grid) = ", -maxi, maxi)
print("jrange (grid) = ", 0, maxj)
print("krange (grid) = ", 0, maxk)
print("height_wg (grid) = ", grid_height_wg)
print("height_bsio2 (grid) = ", grid_height_bsio2)
print("hwidth_wg (grid) = ", grid_hwidth_wg)
print("length_wg1 (grid) = ", grid_length_wg1)
print("length_mmi (grid) = ", grid_length_mmi)
print("hwidth_mmi (grid) = ", grid_hwidth_mmi)
print("hwidth_sep (grid) = ", grid_hwidth_sep)

print("irange (real) = ", -maxi*real_dx, maxi*real_dx)
print("jrange (real) = ", 0, maxj*real_dx)
print("krange (real) = ", 0, maxk*real_dx)


--- create dielectric structure

scene1 = Scene{value=1.}
box_bsio2 = Box{ 
   from={-maxi-1,-1,-1},
   to={maxi+1,grid_height_bsio2,maxk+1}
}

box_wg1 = Box{
   from={-grid_hwidth_wg,grid_height_bsio2,-1},
   to={grid_hwidth_wg,grid_height_bsio2+grid_height_wg,grid_length_wg1 }
}

box_mmi = Box{
   from={-grid_hwidth_mmi,grid_height_bsio2,grid_length_wg1},
   to={grid_hwidth_mmi,grid_height_bsio2+grid_height_wg,grid_length_wg1+grid_length_mmi }
}

i0 = -grid_hwidth_sep-grid_hwidth_wg
i1 = -grid_hwidth_sep+grid_hwidth_wg

box_wg2 = Box{
   from={i0,grid_height_bsio2,grid_length_wg1+grid_length_mmi},
   to={i1,grid_height_bsio2+grid_height_wg,maxk+1 }
}

i0 = grid_hwidth_sep-grid_hwidth_wg
i1 = grid_hwidth_sep+grid_hwidth_wg

box_wg3 = Box{
   from={i0,grid_height_bsio2,grid_length_wg1+grid_length_mmi},
   to={i1,grid_height_bsio2+grid_height_wg,maxk+1 }
}

scene1:add{ box_bsio2, depth=1, value=eps_sio2 }
scene1:add{ box_wg1, depth=1, value=eps_si }
scene1:add{ box_mmi, depth=1, value=eps_si }
scene1:add{ box_wg2, depth=1, value=eps_si }
scene1:add{ box_wg3, depth=1, value=eps_si }
grid1 = Grid{from={-maxi,grid_height_bsio2-10,0},to={maxi,grid_height_bsio2+grid_height_wg+1,maxk}}
gridv =  Grid{yee=false,from={-maxi,grid_height_bsio2-10,0},to={maxi,grid_height_bsio2+grid_height_wg+1,maxk},offset={-maxi,0,0},cells={50,100,100}}

cfg:CREATE_GEO{"scene1", scene=scene1, grid=grid1, method="default",comps=3, silent=false, on=false }
cfg:CREATE_PREVIEW{"scene1", scene=scene1, grid=gridv, method="default", silent=false, on=true }


--- GRID Definition

cfg:GRID{

   dim = 3,
   partition = { 0, 1 },
   ncyc = 1000,
   dt = 0.576,
   irange = { -maxi, maxi },
   jrange = { 0, maxj },
   krange = { 0, maxk }

}
--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 BOX{
	    { -maxi-1, maxi+1, 1, 0, grid_height_bsio2-1, 1, 0, maxk, 1, ":", eps_sio2, eps_sio2, eps_sio2 },
	    { -maxi-1, maxi+1, 1, grid_height_bsio2, maxj, 1, 0, maxk, 1, ":", eps_si, eps_si, eps_si }
	 },
	 LOAD_GEO{ "scene1" }
      },
      on = true
   },

   OUT{
      file = { "VTK", "test_xy_en" },
      type = { "En", "N" },
      time = { 0, 1000, 500 },
      REG{
	 BOX{
	    { -50,50, 1,-50, 50, 1, 0, 0, 1 }
	 }
      }
   },

   OUT{
      file = { "VTK", "test_xy_e" },
      type = { "E", "N" },
      time = { 0, 1000, 500 },
      REG{
	 BOX{
	    { -50,50, 1,-50, 50, 1, 0, 0, 1 }
	 }
      }
   },

   OUT{
      file = { "VTK", "test_xz_e" },
      type = { "E", "N" },
      time = { 0, 1000, 500 },
      REG{
	 BOX{
	    { -50, 50, 1, 0, 0, 1, -50, 50, 1 }
	 }
      }
   },

   OUT{
      file = { "VTK", "test_yz_e" },
      type = { "E", "N" },
      time = { 0, 1000, 500 },
      REG{
	 BOX{
	    { 0, 0, 1, -50, 50, 1, -50, 50, 1 }
	 }
      }
   },

}

--- BOUND Definition

cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },

   PML{
      cells = 11,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

--- SRC Definition(s)

cfg:SRC{
   TFSFBOX{
      invlambda = 0.0333333333,
      amplitude = 1.,
      pulse = { 
	 shape="Gaussian", 
	 width=400,
	 offset=0, attack=800, sustain=400, decay=800 
      },
      planewave = { phi=0, theta=0.0, psi=90.0, nrefr=1.0 }
   },
   REG{
      BOX{
	 { -35, 35, 1, -35, 35, 1, -35, 35, 1 }
      }
   },
   on = false
}


--- CREATE: config.<part>.in
cfg:CREATE()