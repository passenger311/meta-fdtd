cfg:DIAG{
   PSPEC{
      file = "fft_t",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, hdist_tfsf_kp-2, hdist_tfsf_kp-2, 1 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_r",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=0, psi=ppsi }
   },
   REG{
      BOX{
         { imin, imax, step_fft, jmin, jmax, step_fft, -hdist_ntff_kn, -hdist_ntff_kn, 1 }
      }
   },
   on = diagpspec
}

