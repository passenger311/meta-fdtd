include(basic.m4)
include(modules.m4)
include(fields.m4)
include(regloop.m4)

define({M4_CALC_FILTER_COMPONENT},{
M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
  out_value = $1$2(i,j,k) * diag%beta
  do q = 0,diag%p-1,1
    mod_pos = MOD( q + diag%h_pos, diag%p )
    out_value = out_value + &
      (diag%buf_$1(mod_pos,i,j,k,$3)*diag%alpha(q))
  enddo !q
  diag%energy_density = diag%energy_density + &
    ( REALPART(out_value)**2 + IMAGPART(out_value)**2 )
  diag%buf_$1(diag%h_pos,i,j,k,$3) = out_value})
})

define({M4_CALC_FILTER_FIELD},{
M4_CALC_FILTER_COMPONENT($1,x,0)
M4_CALC_FILTER_COMPONENT($1,y,1)
M4_CALC_FILTER_COMPONENT($1,z,2)
})

