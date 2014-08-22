--- DFT

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy_middle",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { imin, imax, step_dft, jmin, jmax, step_dft, hdist_tfsf_kp-2, hdist_tfsf_kp-2, 1 }
      }
   },
   on = false
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xz0",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { imin, imax, step_dft, 0, 0, 1, kmin, kmax, step_dft }
--         { -hdist_ntff_i, hdist_ntff_i, step_dft, 0, 0, 1, -hspacer-2*cyl_diameter, hspacer+2*cyl_diameter , step_dft }
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xz1",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { imin, imax, step_dft, 20, 20, 1, kmin, kmax, step_dft }
--         { -hdist_ntff_i, hdist_ntff_i, step_dft, 0, 0, 1, -hspacer-2*cyl_diameter, hspacer+2*cyl_diameter , step_dft }
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_yz",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { 0, 0, 1, jmin, jmax, step_dft, kmin, kmax, step_dft }
--         { 0, 0, 1, -hdist_ntff_j, hdist_ntff_j, step_dft, -hspacer-2*cyl_diameter, hspacer+2*cyl_diameter , step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         { imin, imax, step_dft, jmin, jmax, step_dft, hspacer, hspacer, step_dft }
--         { 0, 0, 1, -hdist_ntff_j, hdist_ntff_j, step_dft, -hspacer-2*cyl_diameter, hspacer+2*cyl_diameter , step_dft }
      }
   }
}

