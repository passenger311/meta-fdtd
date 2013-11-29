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
      outfile = "F_xy0",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         box_xy0
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy1",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         box_xy1
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         box_xy2
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy3",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         box_xy3
      }
   },
}
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "F_xy4",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_dft },
      mode = "F"
   },
   REG{
      BOX{
         box_xy4
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
      }
   },
   on = false
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
         { imin-size_pml, imax+size_pml, step_dft, jmin-size_pml, jmax+size_pml, step_dft, 0, 0, 1}
      }
   },
   on =false
}
