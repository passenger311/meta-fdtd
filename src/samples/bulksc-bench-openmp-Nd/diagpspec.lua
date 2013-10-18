cfg:DIAG{
   PSPEC{
      file = "fft_T",
      reffile = "fft_ref2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { size_tfsf-size_pad, 0, 0 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_R",
      reffile = "fft_ref2",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { imin, 0, 0 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_t",
      reffile = "fft_ref1",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { size_tfsf-size_pad, 0, 0 }
      }
   },
   on = diagpspec
}
cfg:DIAG{
   PSPEC{
      file = "fft_r",
      reffile = "fft_ref1",
      time = { ncyc_probe_start, ncycles, (ncycles-ncyc_probe_start+1)/samp_fft },
      phasewrap = { 1, 1 },
      mode = "Eap",
      polarize = { phi=0, theta=90, psi=ppsi }
   },
   REG{
      POINT{
         { imin, 0, 0 }
      }
   },
   on = diagpspec
}
