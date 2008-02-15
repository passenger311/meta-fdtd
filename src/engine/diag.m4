include(basic.m4)
include(helper.m4)

define({M4_DIAGLIST}, {esyscmd({./mlist.pl diag.def 2>/dev/null})})
define({M4_FOREACH_DIAG},M4_FOREACHQ({XXX},{M4_DIAGLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_DIAG2},M4_FOREACHQ({XXX},{M4_DIAGLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))