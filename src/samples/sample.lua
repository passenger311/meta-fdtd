
cfg = CONFIG{scenes=false}

--- GRID Definition

cfg:GRID{

   dim = 3,
   partition = { 0, 1 },
   ncyc = 1000,
   dt = 0.376,
   irange = { -50, 50 },
   jrange = { -50, 50 },
   krange = { -50, 50 }

}


--- CREATE SCENE #1

scene1 = Scene{value=0.}
sphere1 = Sphere{at={0,0,0}, radius = 20 }
sphere2 = Sphere{at={0,0,0}, radius = 18 }
shell1 = BinaryAndNot{sphere1,sphere2}
scene1:add{ shell1, depth=1 }
grid1 = Grid{from={-21,-21,-21},to={21,21,21}}

cfg:CREATE_GEO{"scene1", scene=scene1, grid=grid1, method="default", 
		 comps=3, silent=false }

scene2 = Scene{value=1.}
sphere1 = Sphere{at={0,0,0}, radius = 18 }
scene2:add{ sphere1, depth=1, value=2.0 }
grid2 = Grid{from={-21,-21,-21},to={21,21,21}}

cfg:CREATE_GEO{"scene2", scene=scene2, grid=grid2, method="default", 
		 comps=3, silent=false }

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 LOAD_GEO{ "scene2" }
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
   on = true
}

--- MAT Definition(s)

cfg:MAT{
   DRUDE{
      invlambdapl = 0.0471404520791,
      gammapl = 1.e-3,
      order = 2
   },
   REG{
      LOAD_GEO{ "scene1" }
   }
}


--- CREATE: config.<part>.in
cfg:CREATE()