define({M4_OUTGPLLIST}, {esyscmd({mlist.pl outgpl.def 2>/dev/null})})
define({M4_FOREACH_OUTGPL},M4_FOREACHQ({XXX},{M4_OUTGPLLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_OUTGPL2},M4_FOREACHQ({XXX},{M4_OUTGPLLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))