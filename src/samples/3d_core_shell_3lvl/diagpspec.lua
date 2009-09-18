--- FFT

cfg:DIAG{
   PSPEC{
      file = "fft+xabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { hdist_tfsf_i-2, hdist_tfsf_i-2, 1, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft, -hdist_tfsf_k+2, hdist_tfsf_k-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+yabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, hdist_tfsf_j-2, hdist_tfsf_j-2, 1, -hdist_tfsf_k+2, hdist_tfsf_k-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+zabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft, hdist_tfsf_k-2, hdist_tfsf_k-2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-xabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, -hdist_tfsf_i+2, 1, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft, -hdist_tfsf_k+2, hdist_tfsf_k-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-yabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, -hdist_tfsf_j+2, -hdist_tfsf_j+2, 1, -hdist_tfsf_k+2, hdist_tfsf_k-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-zabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft, -hdist_tfsf_k+2, -hdist_tfsf_k+2, 1 }
      }
   }
}


cfg:DIAG{
   PSPEC{
      file = "fft+xscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { hdist_ntff_i+2, hdist_ntff_i+2, 1, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft, -hdist_ntff_k-2, hdist_ntff_k+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+yscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, hdist_ntff_j+2, hdist_ntff_j+2, 1, -hdist_ntff_k-2, hdist_ntff_k+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+zscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft, hdist_ntff_k+2, hdist_ntff_k+2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-xscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, -hdist_ntff_i-2, 1, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft, -hdist_ntff_k-2, hdist_ntff_k+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-yscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=0.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, -hdist_ntff_j-2, -hdist_ntff_j-2, 1, -hdist_ntff_k-2, hdist_ntff_k+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-zscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=0, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft, -hdist_ntff_k-2, -hdist_ntff_k-2, 1 }
      }
   }
}
