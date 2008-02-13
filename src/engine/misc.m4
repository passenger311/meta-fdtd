define({M4_WRITE_DBG}, {M4_IFELSE_DBG({write(6,*) "!DBG (",TRIM(modname),") ", $1})})
define({M4_WRITE_INFO}, {write(STDOUT,*) "!INF (",TRIM(modname),") ", $1 })
define({M4_WRITE_WARN}, {write(STDOUT,*) "!WRN (",TRIM(modname),") ", $1 })
define({M4_FATAL_ERROR}, {write(STDERR,*) "!ERR (",TRIM(modname),") ", $1
stop
})
define({M4_SYNTAX_ERROR},{
if ( $1 ) then
M4_FATAL_ERROR({"SYNTAX ERROR @ LINE: ",$2})
endif
})
define({M4_PARSE_ERROR},{
if ( $1 ) then
M4_FATAL_ERROR({"PARSE ERROR @ LINE: ",$2})
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
define({M4_READ_ERROR}, {if ( $1 .ne. 0 ) then
M4_FATAL_ERROR({"READ ERROR ",$2})
endif
})
define({M4_READ_DECL},{
    character(len=STRLNG) :: m4_skiptill, m4_linestr, m4_string
    integer :: m4_ios
})
define({M4_READ_LINE},{ # $1->unit, $2->line num, $3 ->arg list
$2 = $2 + 1
read($1,*, iostat = m4_ios) $4
if ( m4_ios .ne. 0 ) then
M4_FATAL_ERROR({"PARSE ERROR @LINE = ",TRIM(i2str($2))})
})
define({M4_READ_ENTER},{ # $1->unit, $2->line num, $3->enter token, $4->expected token

    if ( $3 .ne. $4 ) then
       M4_FATAL_ERROR({"BAD SECTION LABEL @LINE = ",TRIM(i2str($2))})
    endif

})
define({M4_READ_LOOP},{


})
