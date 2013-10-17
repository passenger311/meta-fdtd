include(basic.m4)
include(modules.m4)
include(fields.m4)
include(regloop.m4)

define({M4_PK}, { mat%Pk($1,$2,$3,$4) }  ) 
define({M4_QSUM}, { mat%Qsum($1,$2,$3) }  ) 
define({M4_PSUM}, { mat%Psum($1,$2,$3) }  ) 

