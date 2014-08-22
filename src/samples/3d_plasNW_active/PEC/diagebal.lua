cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         { imin, imax, 1, jmin, jmax, 1, kmin, kmax, 1 }
      }
   },
   OUT{
      file = { "GPL", "ebal_l" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
         { imin, imax, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_ldiff" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
         { imin, imax, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      }
   },
      OUT{
      file = { "GPL", "ebal_lds" },
      type = { "DS", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
         { imin, imax, 1, jmin, jmax, 1, kmin, kmax, 1 }
         }
      }
   },

   on = diagebal
}
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         { -40, 40, 1, jmin, jmax, 1, kmin, 40, 1 }
      }
   },
   OUT{
      file = { "GPL", "ebal_s" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { -9, 9, 1, -cyl_length-5, cyl_length+5, 1, kmin, 16, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_sdiff" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            { -9, 9, 1, -cyl_length-5, cyl_length+5, 1, kmin, 16, 1 }
         }
      }
   },
      OUT{
      file = { "GPL", "ebal_sds" },
      type = { "DS", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            { -9, 9, 1, -cyl_length-5, cyl_length+5, 1, kmin, 16, 1 }
         }
      }
   },

   on = diagebal
}
--[[
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
      }
   },
   OUT{
      file = { "GPL", "ebal_xy0" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_xy0_diff" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            { imin, imax, 1, 0, 0, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },
   on = diagebal
}
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      LOAD_GEO{ "xy0" }
   },
   OUT{
      file = { "GPL", "ebal_xy0_en" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         LOAD_GEO{ "xy0" }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff0_en" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            LOAD_GEO{ "xy0" }
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
         { imin, imax, 1, 20, 20, 1, kmin, kmax, 1 }
      }
   },
   OUT{
      file = { "GPL", "ebal_xy1" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 1, 20, 20, 1, kmin, kmax, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_xy1_diff" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            { imin, imax, 1, 20, 20, 1, kmin, kmax, 1 }
         }
      },
      on = true
   },
   on = diagebal
}
cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      LOAD_GEO{ "xy1" }
   },
   OUT{
      file = { "GPL", "ebal_xy1_en" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         LOAD_GEO{ "xy1" }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff1_en" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
            LOAD_GEO{ "xy1" }
         }
      },
      on = false
   },
   on = diagebal
}
--]]
