include(basic.m4)
include(modules.m4)
include(fields.m4)
include(regloop.m4)

define({M4_CALC_FILTER_LOOP},{
out_value = $1$2(i,j,k) * diag%beta
  do q = 0,diag%p-1,1
    mod_pos = MOD( q + diag%h_pos_$1, diag%p )
    out_value = out_value + &
      (diag%buf_$1(mod_pos,i,j,k,$3)*diag%alpha(q))
  enddo !q
  diag%buf_$1(diag%h_pos_$1,i,j,k,$3) = out_value
})

define({M4_CALC_FILTER_COMPONENT},{
M4_REGLOOP_EXPR(diag%reg_outset,p,i,j,k,w,{M4_CALC_FILTER_LOOP($1,$2,$3)})
})

define({M4_CALC_FILTER_FIELD},{
M4_CALC_FILTER_COMPONENT($1,x,0)
M4_CALC_FILTER_COMPONENT($1,y,1)
M4_CALC_FILTER_COMPONENT($1,z,2)
})

define({M4_CALC_FILTER_COMPONENT_OF_MAT},{
M4_REGLOOP_EXPR(mat_reg,p,i,j,k,w,{
if ( ( i .ge. reg%is ) .and. ( i .le. reg%ie ) .and. &
    ( j .ge. reg%js ) .and. ( j .le. reg%je ) .and. &
    ( k .ge. reg%ks ) .and. ( k .le. reg%ke ) ) then
if ( diag%mask(i,j,k) ) then
  out_value = mat%P$1(p,m) * diag%beta
  do q = 0,diag%p-1,1
    mod_pos = MOD( q + diag%h_pos_P, diag%p )
    out_value = out_value + &
      (diag%buf_P(mod_pos,i,j,k,$2,$3)*diag%alpha(q))
  enddo !q
  diag%buf_P(diag%h_pos_P,i,j,k,$2,$3) = out_value
endif
endif})
})

define({M4_CALC_FILTER_FIELD_OF_MAT},{
M4_CALC_FILTER_COMPONENT_OF_MAT(x,0,$1)
M4_CALC_FILTER_COMPONENT_OF_MAT(y,1,$1)
M4_CALC_FILTER_COMPONENT_OF_MAT(z,2,$1)
})


