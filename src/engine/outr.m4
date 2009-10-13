include(basic.m4)
include(version.m4)
include(helper.m4)
include(regloop.m4)

define({M4_OUTRLIST}, {esyscmd({./scripts/mlist.lua outr.def 2>/dev/null})})
define({M4_FOREACH_OUTR},M4_FOREACHQ({XXX},{M4_OUTRLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_OUTR2},M4_FOREACHQ({XXX},{M4_OUTRLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))