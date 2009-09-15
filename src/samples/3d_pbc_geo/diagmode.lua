--- DFT

cfg:DIAG{
   MODE{
      file = "invlambda2.in",
      outfile = "EHT",
      time = { 0, ncycles, (ncycles+1)/samp_dft },
      mode = "EHT"
   },
   REG{
      BOX{
         { -hdist_ntff_i+1, hdist_ntff_i, 1, -hdist_ntff_j+1, hdist_ntff_j, 1, hsio2height+2*hholeheight, hsio2height+2*hholeheight, 1 }
      }
   }
}
