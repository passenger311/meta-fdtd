[ Ex,Ey,Hx,Hy,range] = wgmode3k( 'geo_inj.in' ,0.0072046109510086,1,1,2.7);
wgsave3k( 'tfsfex.set', Ex, Hy, range, range );
wgsave3k( 'tfsfey.set', Ey, -Hx, range, range );
quit;
