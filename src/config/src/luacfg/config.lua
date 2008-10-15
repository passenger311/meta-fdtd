
cfg = CONFIG()

--- GRID Definition

cfg:GRID{

   dim = 2,
   partition = { 0, 1 },
   ncyc = 9000,
   dt = 0.706,
   irange = { -200, 200 },
   jrange = { -200, 200 },
   krange = { 20, 20 }

}


--- CREATE SCENE #1

scene1 = Scene()
scene1:add{ Sphere{at={0,0,0},radius = 100.} }
grid1 = Grid{from={-101,-101,0},to={101,101,0}}

--- FDTD Definition

cfg:FDTD{

   EPSILON{
      REG{
	 GEO { scene = scene1, grid = grid1 } 
      }
   },

   OUT{
      file = { "VTK", "test_xy_en" },
      type = { "En", "N" },
      REG{
	 BOX{
	    { -200, 200, 1,-200, 200, 1 }, {}
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


--- MAT Definition(s)

cfg:MAT{
   
   DRUDE{
      invlambdapl = 0.0471404520791,
      gammapl = 1.e-10,
      order = 1
   },
   REG{
      BOX{
	 { -30, 30, 1, -30, 30, 1 }, { 1.,1.,1. },
	 { -30, -30, 1, -30, 30, 1 }, { 1.,.5,.5 },
	 { 30, 30, 1, -30, 30, 1 }, { 0.,.5,.5 },
	 { -30, 30, 1, -30, -30, 1 }, { .5,1.,.5 },
	 { -30, 30, 1, 30, 30, 1 }, { .5,0.,.5 },
	 { -30, -30, 1, -30, -30, 1 }, { .5,.5,.25 },
	 { -30, -30, 1, 30, 30, 1 }, { .5,0.,.25 },
	 { 30, 30, 1, -30, -30, 1 }, { 0.,.5,.25 },
	 { 30, 30, 1, 30, 30, 1 }, { 0.,0.,.25 },
      }
   },
   on = true
}

cfg:MAT{
   
   BLOCH{
      invlambdal = 0.0471404520791,
      gammal = 1.e-10,
      dipole = { {2.,0.},{2.,0.},{1,0.} },
      carrier = { 1. , 0.5 },
      gammanr = 1.e-9,
      pump = 1.,
      satmodel = 0
   },
   REG{
      BOX{
	 { -30, 30, 1, -30, 30, 1 }, { 1.,1.,1. },
	 { -30, -30, 1, -30, 30, 1 }, { 1.,.5,.5 },
	 { 30, 30, 1, -30, 30, 1 }, { 0.,.5,.5 },
	 { -30, 30, 1, -30, -30, 1 }, { .5,1.,.5 },
	 { -30, 30, 1, 30, 30, 1 }, { .5,0.,.5 },
	 { -30, -30, 1, -30, -30, 1 }, { .5,.5,.25 },
	 { -30, -30, 1, 30, 30, 1 }, { .5,0.,.25 },
	 { 30, 30, 1, -30, -30, 1 }, { 0.,.5,.25 },
	 { 30, 30, 1, 30, 30, 1 }, { 0.,0.,.25 },
      }
   },
   on = true
}

--- CREATE: config.<part>.in
cfg:CREATE()