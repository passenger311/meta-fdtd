--[[
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         box_xy0
      }
   },
   OUT{
      file = { "GPL", "ebal_xy1" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy0
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff1" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 4 } ,
      REG{
         BOX{
            box_xy0
         }
      },
      on = false
   },
   on = diagebal
}
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         box_xy3
      }
   },
   OUT{
      file = { "GPL", "ebal_xy3" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy3
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff3" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 4 } ,
      REG{
         BOX{
            box_xy3
         }
      },
      on = false
   },
   on = diagebal
}
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         box_xy4
      }
   },
   OUT{
      file = { "GPL", "ebal_xy4" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy4
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff4" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 4 } ,
      REG{
         BOX{
            box_xy4
         }
      },
      on = false
   },
   on = diagebal
}
--]]
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         box_xy2
      }
   },
   OUT{
      file = { "GPL", "ebal_xy2" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         BOX{
            box_xy2
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff2" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 4 } ,
      REG{
         BOX{
            box_xy2
         }
      },
      on = false
   },
   on = diagebal
}

cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      LOAD_GEO{ "xy2" }
   },
   OUT{
      file = { "GPL", "ebal_xy2_en" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 4 },
      REG{
         LOAD_GEO{ "xy2" }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff2_en" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 4 } ,
      REG{
         BOX{
            LOAD_GEO{ "xy2" }
         }
      },
      on = false
   },
   on = diagebal
}
