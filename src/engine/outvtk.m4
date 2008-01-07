define({M4_OUTVTKLIST}, {esyscmd({mlist.pl outgpl.def 2>/dev/null})})
define({M4_FOREACH_OUTVTK},M4_FOREACHQ({XXX},{M4_OUTVTKLIST},{format({%s%s%s},$1,{XXX},$2)}))
define({M4_FOREACH_OUTVTK2},M4_FOREACHQ({XXX},{M4_OUTVTKLIST},{format({%s%s%s%s%s},$1,{XXX},$2,{XXX},$3)}))