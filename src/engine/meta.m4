changequote(`{', `}')

include({helper.m4})
include({reglist.m4})
include({mat.m4})
include({diag.m4})
include({outgpl.m4})
include({modules.m4})
include({dim.m4})
include({fdtd.m4})
include({misc.m4})

define({M4_COMMENT},{})

define({M4_FTYPE}, {ifdef({M4_CF}, {complex(kind=8)}, {real(kind=8)})})

define({M4_MPICOORD},{ifelse("M4_SDIM","3",{IMIN,JMIN,K$1},{ifelse("M4_DIM","2",{IMIN,J$1,0},{I$1,0,0})})})
  
define({M4_IS3D},{ifelse("M4_SDIM","3", {.true.}, {.false.})})
define({M4_IS2D},{ifelse("M4_SDIM","2", {.true.}, {.false.})})
define({M4_IS1D},{ifelse("M4_SDIM","1", {.true.}, {.false.})})
define({M4_ISDBG},{ifdef({M4_DBG}, {.true.}, {.false.})})
define({M4_ISMPI}, {ifdef({M4_MPI}, {.true.}, {.false.})})
define({M4_ISMPELOG}, {ifdef({M4_MPELOG}, {.true.}, {.false.})})
define({M4_ISOMP}, {ifdef({M4_OMP}, {.true.}, {.false.})})
define({M4_ISCF}, {ifdef({M4_CF}, {.true.}, {.false.})})
define({M4_ISWMU},{ifdef({M4_WMU},{.true.},{.false.})})
define({M4_ISTE}, {ifdef({M4_TE}, {.true.}, {.false.})})
define({M4_ISNG}, {ifdef({M4_NG}, {.true.}, {.false.})})

define({M4_IFELSE_DBG}, {ifdef({M4_DBG}, {$1}, {$2})})
define({M4_IFELSE_MPI}, {ifdef({M4_MPI}, {$1}, {$2})})
define({M4_IFELSE_MPELOG}, {ifdef({M4_MPELOG}, {$1}, {$2})})
define({M4_IFELSE_OMP}, {ifdef({M4_OMP}, {$1}, {$2})})
define({M4_IFELSE_CF}, {ifdef({M4_CF}, {$1}, {$2})})
define({M4_IFELSE_WMU},{ifdef({M4_WMU},{$1},{$2})})
define({M4_IFELSE_TE}, {ifdef({M4_TE}, {$1}, {$2})})
define({M4_IFELSE_NG}, {ifdef({M4_NG}, {$1}, {$2})})
define({M4_IFELSE_ARCH},{ifelse("M4_ARCH","$1",$2,$3)}) 

define({M4_MPE_SECTION}, { M4_IFELSE_MPELOG({
call StartMPELog($1)  
$2
call StopMPELog($1)  
},$2)})
