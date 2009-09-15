--- FFT

cfg:DIAG{
   PSPEC{
      file = "fft+zabs",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+1, hdist_tfsf_i, step_fft, -hdist_tfsf_j+1, hdist_tfsf_j, step_fft, hdist_tfsf_k-2, hdist_tfsf_k-2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-zabs",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+1, hdist_tfsf_i, step_fft, -hdist_tfsf_j+1, hdist_tfsf_j, step_fft, -hdist_tfsf_k+2, -hdist_tfsf_k+2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft+zscat",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i+1, hdist_ntff_i, step_fft, -hdist_ntff_j+1, hdist_ntff_j, step_fft, hdist_ntff_k+2, hdist_ntff_k+2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-zscat",
      time = { 0, ncycles, (ncycles+1)/1024 },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i+1, hdist_ntff_i, step_fft, -hdist_ntff_j+1, hdist_ntff_j, step_fft, -hdist_ntff_k-2, -hdist_ntff_k-2, 1 }
      }
   }
}
