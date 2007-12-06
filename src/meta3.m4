include(`reglist.m4')

define(`M4_FTYPE', `ifdef(`M4_CMPLX', `complex(kind=8)', `real(kind=8)')')
define(`M4_FTYPE_NUM', `ifdef(`M4_CMPLX', `2', `1')')
