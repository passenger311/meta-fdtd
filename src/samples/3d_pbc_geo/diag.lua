cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "EHT",
      time = { 0, ncycles, (ncycles+1)/2048 },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, 5, -hdist_ntff_j, hdist_ntff_j, 5, hsio2height+2*hholeheight, hsio2height+2*hholeheight, 1 }
      }
   }
}

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
         { -hdist_tfsf_i, hdist_tfsf_i, step_fft, -hdist_tfsf_j, hdist_tfsf_j, step_fft, hdist_tfsf_k-2, hdist_tfsf_k-2, 1 }
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
         { -hdist_tfsf_i, hdist_tfsf_i, step_fft, -hdist_tfsf_j, hdist_tfsf_j, step_fft, -hdist_tfsf_k+2, -hdist_tfsf_k+2, 1 }
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
         { -hdist_ntff_i, hdist_ntff_i, step_fft, -hdist_ntff_j, hdist_ntff_j, step_fft, hdist_ntff_k+2, hdist_ntff_k+2, 1 }
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
         { -hdist_ntff_i, hdist_ntff_i, step_fft, -hdist_ntff_j, hdist_ntff_j, step_fft, -hdist_ntff_k-2, -hdist_ntff_k-2, 1 }
      }
   }
}
