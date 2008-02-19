changequote(`{', `}')

define({M4_COMMENT},{})
define({M4_FTYPE}, {ifdef({M4_CF}, {complex(kind=8)}, {real(kind=8)})})
define({M4_SDIM},{ifelse("M4_DIM","3",{3$1},ifelse("M4_DIM","2",{2$1},{1$1}))})
define({M4_MPICOORD},{ifelse("M4_SDIM","3",{IMIN,JMIN,K$1},{ifelse("M4_DIM","2",{IMIN,J$1,0},{I$1,0,0})})})
define({M4_IS3D},{ifelse("M4_SDIM","3", {.true.}, {.false.})})
define({M4_IS2D},{ifelse("M4_SDIM","2", {.true.}, {.false.})})
define({M4_IS1D},{ifelse("M4_SDIM","1", {.true.}, {.false.})})
define({M4_ISDBG},{ifdef({M4_DBG}, {.true.}, {.false.})})
define({M4_ISMPI}, {ifdef({M4_MPI}, {.true.}, {.false.})})
define({M4_ISMPELOG}, {ifdef({M4_MPELOG}, {.true.}, {.false.})})
define({M4_ISOMP}, {ifdef({M4_OMP}, {.true.}, {.false.})})
define({M4_ISCF}, {ifdef({M4_CF}, {.true.}, {.false.})})
define({M4_CONJ},{ifdef({M4_CF}, {conjg($1)}, {$1})})
define({M4_ISWMU},{ifdef({M4_WMU},{.true.},{.false.})})
define({M4_ISTE}, {ifdef({M4_TE}, {.true.}, {.false.})})
define({M4_ISNG}, {ifdef({M4_NG}, {.true.}, {.false.})})
define({M4_IFELSE_DBG}, {ifdef({M4_DBG}, {$1}, {$2})})
define({M4_IFELSE_HD5}, {ifdef({M4_HD5}, {$1}, {$2})})
define({M4_IFELSE_MPI}, {ifdef({M4_MPI}, {$1}, {$2})})
define({M4_IFELSE_MPELOG}, {ifdef({M4_MPELOG}, {$1}, {$2})})
define({M4_IFELSE_OMP}, {ifdef({M4_OMP}, {$1}, {$2})})
define({M4_IFELSE_CF}, {ifdef({M4_CF}, {$1}, {$2})})
define({M4_IFELSE_WMU},{ifdef({M4_WMU},{$1},{$2})})
define({M4_IFELSE_TEPS}, {ifdef({M4_TEPS}, {$1}, {$2})})
define({M4_IFELSE_NG}, {ifdef({M4_NG}, {$1}, {$2})})
define({M4_IFELSE_ARCH},{ifelse("M4_ARCH","$1",$2,$3)}) 
define({M4_IFELSE_3D},{ifelse("M4_SDIM","3", {$1}, {$2})})
define({M4_IFELSE_2D},{ifelse("M4_SDIM","2", {$1}, {$2})})
define({M4_IFELSE_1D},{ifelse("M4_SDIM","1", {$1}, {$2})})
define({M4_IFELSE_TE}, {M4_IFELSE_2D({ifdef({M4_TE}, {$1}, {$2})},{$1})})
define({M4_IFELSE_TM}, {M4_IFELSE_2D({ifdef({M4_TM}, {$1}, {$2})},{$1})})
define({M4_DIM123},{M4_IFELSE_1D({$1},{M4_IFELSE_2D({$2},{$3})})})

define({M4_MPE_SECTION}, { M4_IFELSE_MPELOG({
call StartMPELog($1)  
$2
call StopMPELog($1)  
},$2)})

define({M4_WRITE_DBG}, {M4_IFELSE_DBG({write(6,*) "!DBG (",TRIM(modname),") ", $1})})
define({M4_WRITE_INFO}, {write(STDOUT,*) "!INF (",TRIM(modname),") ", $1 })
define({M4_WRITE_WARN}, {write(STDOUT,*) "!WRN (",TRIM(modname),") ", $1 })
define({M4_FATAL_ERROR}, {write(STDERR,*) "!ERR (",TRIM(modname),") ", $1
stop
})
define({M4_SYNTAX_ERROR},{
if ( $1 ) then
M4_FATAL_ERROR({"INPUT ERROR (@",TRIM(i2str($2-1)),") EXPECTED ",$3})
endif
})
define({M4_PARSE_ERROR},{
if ( $1 ) then
M4_FATAL_ERROR({"INPUT ERROR (@",TRIM(i2str($2-1)),") $3"})
endif
})
define({M4_ALLOC_ERROR}, {if ( $1 .ne. 0 ) then
M4_FATAL_ERROR({"OUT OF MEMORY ",$2})
endif
})
define({M4_OPEN_ERROR}, {if ( $1 .ne. 0 ) then
M4_FATAL_ERROR({"COULD NOT OPEN FILE ",$2})
endif
})
