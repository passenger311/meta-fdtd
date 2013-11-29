cfg:DIAG{
   PSPEC{
      file = "fft_T",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, hdist_tfsf_kp-2, hdist_tfsf_kp-2, 1 }
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_R",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, -hdist_tfsf_kn, -hdist_tfsf_kn, 1 }
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_R2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, -2*cyl_diameter-hspacer-5, -2*cyl_diameter-hspacer-5, 1 }
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_t",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, hdist_tfsf_kp-2, hdist_tfsf_kp-2, 1 }
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_r",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, -hdist_tfsf_kn, -hdist_tfsf_kn, 1 }
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_r2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, -2*cyl_diameter-hspacer-5, -2*cyl_diameter-hspacer-5, 1, kmin, kmax, step_fft }
      }
   },
   on = false,
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer+cyl_radius, hspacer+cyl_radius, 1, -cyl_length, cyl_length, step_fft}
      }
   },
   on = false
--   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_xy0",
      reffile = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer, hspacer+2*cyl_radius, step_fft, inj_pos-2, inj_pos-2, 1}
      }
   },
--   on = false
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_xy1",
      reffile = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer, hspacer+2*cyl_radius, step_fft, math.floor(inj_pos/2), math.floor(inj_pos/2), 1}
      }
   },
--   on = false
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_xy2",
      reffile = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer, hspacer+2*cyl_radius, step_fft, -2, -2, 1}
      }
   },
--   on = false
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_xy3",
      reffile = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer, hspacer+2*cyl_radius, step_fft, -math.floor(inj_pos/2), -math.floor(inj_pos/2), 1}
      }
   },
--   on = false
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_xy4",
      reffile = "fft_ref",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { -cyl_radius, cyl_radius, step_fft, hspacer, hspacer+2*cyl_radius, step_fft, -inj_pos+2, -inj_pos+2, 1}
      }
   },
--   on = false
   on = diagpspec
}


