
cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "EHN",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHN"
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, 10, -hdist_tfsf_j+2, hdist_tfsf_j-2, 5 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft+xabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { hdist_tfsf_i-2, hdist_tfsf_i-2, step_fft, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+yabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, hdist_tfsf_j-2, hdist_tfsf_j-2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-xabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, -hdist_tfsf_i+2, 1, -hdist_tfsf_j+2, hdist_tfsf_j-2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-yabs",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, step_fft, -hdist_tfsf_j+2, -hdist_tfsf_j+2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft+xscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { hdist_ntff_i+2, hdist_ntff_i+2, 1, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft+yscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, hdist_ntff_j+2, hdist_ntff_j+2, 1 }
      }
   }
}

cfg:DIAG{
   PSPEC{
      file = "fft-xscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=0, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, -hdist_ntff_i-2, 1, -hdist_ntff_j-2, hdist_ntff_j+2, step_fft }
      }
   }
}
cfg:DIAG{
   PSPEC{
      file = "fft-yscat",
      time = { 0, ncycles, (ncycles+1)/sampl_fft },
      phasewrap = { 1, 0 },
      mode = "S",
      polarize = { phi=90, theta=90, psi=90.0 }
   },
   REG{
      BOX{
         { -hdist_ntff_i-2, hdist_ntff_i+2, step_fft, -hdist_ntff_j-2, -hdist_ntff_j-2, 1 }
      }
   }
}

--- DFT

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft+x",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { hdist_ntff_i, hdist_ntff_i, 1, -hdist_ntff_j, hdist_ntff_j, step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft+y",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, hdist_ntff_j, hdist_ntff_j, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft-x",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, -hdist_ntff_i, 1, -hdist_ntff_j, hdist_ntff_j, step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft-y",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, -hdist_ntff_j, -hdist_ntff_j, 1 }
      }
   }

}

