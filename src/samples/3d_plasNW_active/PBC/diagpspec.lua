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
         { imin, imax, step_fft, jmin, jmax, step_fft, -hdist_ntff_kn, -hdist_ntff_kn, 1 }
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
         { imin, imax, step_fft, jmin, jmax, step_fft, -hdist_ntff_kn, -hdist_ntff_kn, 1 }
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
         { imin, imax, step_fft, jmin, jmax, step_fft, -2*cyl_diameter-hspacer-5, -2*cyl_diameter-hspacer-5, 1 }
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
         { -math.ceil(cyl_radius), math.ceil(cyl_radius), step_fft, -cyl_length, cyl_length, 3*step_fft, hspacer+math.ceil(cyl_radius), hspacer+math.ceil(cyl_radius), 1}
      }
   },
--   on = false
   on = diagpspec
}


