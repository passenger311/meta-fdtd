cfg = CONFIG{scenes=true}

cfg:GRID{

   dim = 3,
   partition = { 0, 1 },
   ncyc = 2,
   dt = 0.57,
   irange = { -50, 50 },
   jrange = { -50, 50 },
   krange = { -50, 50 }

}

scene1 = Scene{
   value = 2.0
}
cyl1 = Cylinder{
   at = {0,0,0},
   radius = {30},
   height = {7}
}
sph1 = Sphere{
   at = {0,0,0},
   radius = {30}
}
box1 = Box{
--   at = {-15,0,0},
--   size = {30,30,30}
   from = {-30,-15,-15},
   to = {0,15,15}
}
scene1:add{
   box1,
   value = 3
}
box2 = Box{
--   at = {15,0,0},
--   size = {30,30,30}
   from = {0,-15,-15},
   to = {30,15,15}
}
scene1:add{
   box2,
   value = 4
}

grid1 = Grid{
   from = {-40,-40,-40},
   to = {40,40,40}
}
grid_prev1 = Grid{
   yee = false,
   from = {-40,-40,-40},
   to = {40,40,40}
}

cfg:CREATE_GEO{
   "scene1",
   scene=scene1,
   grid=grid1,
   method="default",
   comps=3,
   silent=false,
   on=true
}

cfg:CREATE_PREVIEW{
   "scene1",
   scene=scene1,
   grid=grid_prev1,
   method="default",
   silent=false,
   on=true
}

cfg:FDTD{
   EPSILON{
      REG{
         LOAD_GEO{ "scene1" }
      }
   },
   OUT{
      file = { "VTK", "Epsx" },
      type = { "Epsx", "N" },
      time = { 0, 0, 1 },
      REG{
         BOX{
            { -40, 40, 1, -40, 40, 1, 0, 0, 1 }
         }
      }
   },
   OUT{
      file = { "VTK", "Epsy" },
      type = { "Epsy", "N" },
      time = { 0, 0, 1 },
      REG{
         BOX{
            { -40, 40, 1, -40, 40, 1, 0, 0, 1 }
         }
      }
   },
   OUT{
      file = { "VTK", "Epsz" },
      type = { "Epsz", "N" },
      time = { 0, 0, 1 },
      REG{
         BOX{
            { -40, 40, 1, -40, 40, 1, 0, 0, 1 }
         }
      }
   }
}

cfg:BOUND{

   config = { 1, 1, 1, 1, 1, 1 },

   PML{
      cells = 10,
      pot = 3.2,
      sigma = 1.94444444444444,
      kappa = 1.1,
   }
}

cfg:SRC{
   TFSFBOX{
      invlambda = 1/240,
      amplitude = 1.0,
      pulse = {
         shape="Gaussian",
         width=200,
         offset=0,
         attack=1000,
         sustain=0,
         decay=2000
      },
      planewave = { phi=90.0, theta=90.0, psi=0.0, nrefr=1 },
      config = {1,1,1,1,1,1}
   },
   REG{
      BOX{
         {-40,40,1,-40,40,1,-40,40,1}
      }
   },
   on = true
}

cfg:CREATE()

