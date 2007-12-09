changequote(`{', `}')

include({helper.m4})
include({reglist.m4})
include({mat.m4})
include({diag.m4})
include({outgpl.m4})


define({M4_ISDBG}, {ifdef({M4_DBG}, {.true.}, {.false.})})
define({M4_ISMPI}, {ifdef({M4_MPI}, {.true.}, {.false.})})
define({M4_ISOMP}, {ifdef({M4_OMP}, {.true.}, {.false.})})
define({M4_ISCF}, {ifdef({M4_CF}, {.true.}, {.false.})})
define({M4_ISTE}, {ifdef({M4_TE}, {.true.}, {.false.})})
define({M4_ISNG}, {ifdef({M4_NG}, {.true.}, {.false.})})

define({M4_IFELSE_DBG}, {ifdef({M4_DBG}, {$1}, {$2})})
define({M4_IFELSE_MPI}, {ifdef({M4_MPI}, {$1}, {$2})})
define({M4_IFELSE_OMP}, {ifdef({M4_OMP}, {$1}, {$2})})
define({M4_IFELSE_CF}, {ifdef({M4_CF}, {$1}, {$2})})
define({M4_IFELSE_TE}, {ifdef({M4_TE}, {$1}, {$2})})
define({M4_IFELSE_NG}, {ifdef({M4_NG}, {$1}, {$2})})

define({M4_ISMPELOG}, {ifdef({M4_MPI}, {ifdef({M4_DBG},{.true.},{})}, {.false.})})
define({M4_IFELSE_MPELOG}, {ifdef({M4_MPI}, {ifdef({M4_DBG},{$1},{})}, {$2})})

define({M4_WRITE_DBG}, {M4_IFELSE_DBG({write(6,*) "!DBG ", $1},{})})
define({M4_FATAL_ERROR}, {write(STDERR,*) "!ERR ", $1
stop
})

define({M4_FTYPE}, {ifdef({M4_CF}, {complex(kind=8)}, {real(kind=8)})})

