include(`helper.m4')
include(`reglist.m4')
include(`mat.m4')
include(`diag.m4')
include(`outgpl.m4')

define(`M4_FTYPE', `ifdef(`M4_CF', `complex(kind=8)', `real(kind=8)')')
