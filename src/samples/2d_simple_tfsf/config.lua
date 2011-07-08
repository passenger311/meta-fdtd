cfg = CONFIG{scenes=true}
-- taskid = select(1,...)

dim = 2
n = { bg = 1, obj = 2 }

-- dofile("scale.lua")
real_wavelength = 1.55 -- [um]
dx = real_wavelength / 25 / n.obj

dt = 0.99* 1/math.sqrt(dim)  
ncycles = 2000
pmls = 11
lambda = real_wavelength / dx

ib = { -100, 100 }
jb = { -100, 100 }

ir = { ib[1]-pmls, ib[2]+pmls }  -- total computational domain
jr = { jb[1]-pmls, jb[2]+pmls }

ir0 = { ir[1]-1, ir[2]+1 }
jr0 = { jr[1]-1, jr[2]+1 }

print( "dx = ", dx )
print( "ir0 = ", ir0[1], ir0[2] )
print( "jr0 = ", jr0[1], jr0[2] )

--- CREATE GEOMETRIES

scene = Scene{ value = n.bg^2 }

sphere = Sphere{ at={0,0,0}, radius = 50.1 }
--scene:add{ sphere, depth = 1, value = n.obj^2 }

grid = Grid{
   from={ir0[1], jr0[1] }, to={ ir0[2], jr0[2]}
}

grid_prev = Grid{
   yee = false,	
   from={ir0[1], jr0[1] }, to={ ir0[2], jr0[2] }
}

cfg:CREATE_GEO{
    "geo",
     scene=scene,
     grid=grid,
}

cfg:CREATE_PREVIEW{
     "geo",
     scene=scene,
     grid=grid_prev
}


--- ASSEMBLE CONFIG

cfg:GRID{
   dim = dim,                       -- number of dimensions
   ncyc = ncycles,                  -- number of time steps
   dt = dt,                         -- time step in units of cell size
   irange = ir,                     -- i range 
   jrange = jr,                     -- j range  
   krange = kr,
--   dx = { conv*1e-9, 1, 1, 1 }
}

cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },  -- 0: pec, 1: pml, 2: pmc, 3: pbc

   PML{
      cells = pmls,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

cfg:FDTD{

   EPSILON{
      REG{
         BOX{
	    { ir0[1], ir0[2], 1,  jr0[1], jr0[2], 1, 
	      ":", n.bg^2, n.bg^2, n.bg^2 } 
         },
         LOAD_GEO{ "geo" }
      },
      on = true
   },

   OUT{
      file = { "VTK", "Ez" },
      type = { "Ez", "N" },
      time = { 0, ncycles, ncycles/100 },
      REG{
         BOX{
            { ir[1], ir[2], 1, jr[1], jr[2], 1 }
         }
      },
      on = true
   },

}


cfg:SRC{
   TFSFBOX{
      invlambda = 1/lambda,
      amplitude = 1.0,
      pulse = { 
         shape="Gaussian",
         width=200,
         offset=0,
         attack=400,
         sustain=2000,
         decay=400
      },
      planewave = { phi=20, theta=90, psi=90, nrefr=n.bg }
   },
   REG{
      BOX{
         { ib[1]+30, ib[2]-30, 1, jb[1]+30, jb[2]-30, 1 }
      }
   },
   on = true
}


--- CREATE: config.<part>.in
cfg:CREATE()
