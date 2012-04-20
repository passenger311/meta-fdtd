--- DFT

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy_middle",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, -hdist_ntff_j, hdist_ntff_j, step_dft, hcoating+hmetal+math.floor(hspacer/2+.5), hcoating+hmetal+math.floor(hspacer/2+.5) , step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy_hole_front",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { -xperf-2, xperf+2, step_dft, -yperf-2, yperf+2, step_dft, hcoating+math.floor(hmetal/2+0.5), hcoating+math.floor(hmetal/2+0.5) , step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy_hole_back",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { -xperf-2, xperf+2, step_dft, -yperf-2, yperf+2, step_dft, hcoating+hspacer+math.floor(3*hmetal/2+0.5), hcoating+hspacer+math.floor(3*hmetal/2+0.5) , step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xz",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, 0, 0, step_dft, -hspacer, 2*hcoating+2*hmetal+2*hspacer , step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_yz",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { 0, 0, 1, -hdist_ntff_j, hdist_ntff_j, 1, -hspacer, 2*hcoating+2*hmetal+2*hspacer , 1 }
      }
   }
}

