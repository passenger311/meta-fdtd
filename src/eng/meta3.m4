include(`helper.m4')
include(`reglist.m4')
include(`mat.m4')
include(`diag.m4')
include(`outgpl.m4')

# --- macro substitutions rules follow

define(`M4_ISDBG', `ifdef(`DBG', `.true.', `.false.')')
define(`M4_ISMPI', `ifdef(`MPI', `.true.', `.false.')')
define(`M4_ISCF', `ifdef(`CF', `.true.', `.false.')')
define(`M4_ISMPI', `ifdef(`MPI', `.true.', `.false.')')
define(`M4_ISTE', `ifdef(`TE', `.true.', `.false.')')
define(`M4_ISNG', `ifdef(`NG', `.true.', `.false.')')

define(`M4_IFELSE_DBG', `ifdef(`DBG', `$1', `$2')')
define(`M4_IFELSE_MPI', `ifdef(`MPI', `$1', `$2')')
define(`M4_IFELSE_CF', `ifdef(`CF', `$1', `$2')')
define(`M4_IFELSE_MPI', `ifdef(`MPI', `$1', `$2')')
define(`M4_IFELSE_TE', `ifdef(`TE', `$1', `$2')')
define(`M4_IFELSE_NG', `ifdef(`NG', `$1', `$2')')

define(`M4_FTYPE', `ifdef(`CF', `complex(kind=8)', `real(kind=8)')')

