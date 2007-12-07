define(`M4_MATLIST', `esyscmd(`cat mat.def 2>/dev/null')')
define(`M4_FOREACH_MAT',M4_FOREACHQ(`XXX',`M4_MATLIST',`format(`%s%s%s',$1,`XXX',$2)'))
define(`M4_FOREACH_MAT2',M4_FOREACHQ(`XXX',`M4_MATLIST',`format(`%s%s%s%s%s',$1,`XXX',$2,`XXX',$3)'))