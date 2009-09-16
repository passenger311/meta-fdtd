--- DFT

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "EHN",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHN"
   },
   REG{
      BOX{
         { -hdist_tfsf_i+2, hdist_tfsf_i-2, 1, 0, 0, 1, -hdist_tfsf_k+2, hdist_tfsf_k-2, 1 }
      }
   }
}

cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft+x",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { hdist_ntff_i, hdist_ntff_i, 1, -hdist_ntff_j, hdist_ntff_j, step_dft, -hdist_ntff_k, hdist_ntff_k, step_dft }
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
         { -hdist_ntff_i, hdist_ntff_i, step_dft, hdist_ntff_j, hdist_ntff_j, 1, -hdist_ntff_k, hdist_ntff_k, step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft+z",
      time = { 0, ncycles, (ncycles+1)/sampl_dft},
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, -hdist_ntff_j, hdist_ntff_j, step_dft, hdist_ntff_k, hdist_ntff_k, 1 }
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
         { -hdist_ntff_i, -hdist_ntff_i, 1, -hdist_ntff_j, hdist_ntff_j, step_dft, -hdist_ntff_k, hdist_ntff_k, step_dft }
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
         { -hdist_ntff_i, hdist_ntff_i, step_dft, -hdist_ntff_j, -hdist_ntff_j, 1, -hdist_ntff_k, hdist_ntff_k, step_dft }
      }
   }
}
cfg:DIAG{
   MODE{
      file = "invlambda.in",
      outfile = "dft-z",
      time = { 0, ncycles, (ncycles+1)/sampl_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i, hdist_ntff_i, step_dft, -hdist_ntff_j, hdist_ntff_j, step_dft, -hdist_ntff_k, -hdist_ntff_k, 1 }
      }
   }
}

