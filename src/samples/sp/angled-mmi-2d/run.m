[ Ex,Ey,Hx,Hy,range] = wgmode3d_2( 'geo_inj.in',0.014409221902017,1,1,2.532);
savetfsfinj_2( 'tfsfex.set', Ex, Hy, range, range );
savetfsfinj_2( 'tfsfey.set', Ey, Hx, range, range );
quit;
