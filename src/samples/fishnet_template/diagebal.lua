cfg:DIAG{
   EBAL{
      time = { 0, ncycles, 1 },
   },
   REG{
      BOX{
         { imin, imax, 1, jmin, jmax, 1, 0, 2*(hcoating+hmetal)+hspacer, 1 }
      }
   },
   OUT{
      file = { "GPL", "ebal" },
      type = { "EnI", "N", ".F." },
      time = { 0, ncycles, 10 },
      REG{
         BOX{
            { imin, imax, 1, jmin, jmax, 1, 0, 2*(hcoating+hmetal)+hspacer, 1 }
         }
      }
   },
   OUT{
      file = { "GPL", "ebal_diff" },
      type = { "En", "N", ".F." },
      time = { 0, ncycles, 10 } ,
      REG{
         BOX{
         { imin, imax, 1, jmin, jmax, 1, 0, 2*(hcoating+hmetal)+hspacer, 1 }
         }
      }
   },
   on = diagebal
}

