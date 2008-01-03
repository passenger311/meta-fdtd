define({M4_SDIM},{ifelse("M4_DIM","3",{3$1},ifelse("M4_DIM","2",{2$1},{1$1}))})
define({M4_IFELSE_3D},{ifelse("M4_SDIM","3", {$1}, {$2})})
define({M4_IFELSE_2D},{ifelse("M4_SDIM","2", {$1}, {$2})})
define({M4_IFELSE_1D},{ifelse("M4_SDIM","1", {$1}, {$2})})
define({M4_DIM123},{M4_IFELSE_1D({$1},{M4_IFELSE_2D({$2},{$3})})})

define({M4_READCOORD},{ifelse("M4_SDIM","3",{$1,$2,$3},{ifelse("M4_DIM","2",{$1,$2},{$1})})})
define({M4_COORD},{ifelse("M4_SDIM","3",{$1,$2,$3},{ifelse("M4_DIM","2",{$1,$2,0},{$1,0,0})})})
define({M4_K},{ifelse("M4_SDIM","3",{$1},{0})})
define({M4_J},{ifelse("M4_SDIM","1",{0},{$1})})
define({M4_RANGE},{ifelse("M4_SDIM","3",{$1,$2,$3},{ifelse("M4_DIM","2",{$1,$2,0:0},{$1,0:0,0:0})})})
define({M4_KRANGE},{ifelse("M4_SDIM","3",{$1},{0:0})})
define({M4_JRANGE},{ifelse("M4_SDIM","1",{0:0},{$1})})
define({M4_KJILOOP},{
M4_IFELSE_3D({		
  do k = $1,$2,$3
    do j = $4,$5,$6 
      do i = $7,$8,$9
},{
  do k = 0, 0, 1
  M4_IFELSE_2D({
  do j = $4,$5,$6
  },{
do j = 0, 0, 1
})
do i = $7,$8,$9
})			
$10
  
enddo
enddo
enddo
})