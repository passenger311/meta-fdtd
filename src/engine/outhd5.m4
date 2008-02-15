include(basic.m4)
include(helper.m4)

define({M4_OUTHD5LIST}, {esyscmd({mlist.pl outgpl.def 2>/dev/null})})
define({M4_FOREACH_OUTHD5},M4_FOREACHQ({XXX},{M4_OUTHD5LIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_OUTHD52},M4_FOREACHQ({XXX},{M4_OUTHD5LIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))