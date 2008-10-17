include(basic.m4)
include(version.m4)
include(helper.m4)
include(regloop.m4)

define({M4_OUTGPLLIST}, {esyscmd({./scripts/mlist.lua outgpl.def 2>/dev/null})})
define({M4_FOREACH_OUTGPL},M4_FOREACHQ({XXX},{M4_OUTGPLLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_OUTGPL2},M4_FOREACHQ({XXX},{M4_OUTGPLLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))