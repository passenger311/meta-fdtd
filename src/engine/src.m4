include(basic.m4)
include(helper.m4)

define({M4_SRCLIST}, {esyscmd({./scripts/mlist.lua src.def 2>/dev/null})})
define({M4_FOREACH_SRC},M4_FOREACHQ({XXX},{M4_SRCLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_SRC2},M4_FOREACHQ({XXX},{M4_SRCLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))