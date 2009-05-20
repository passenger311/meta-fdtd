[ Ex,Ey,Hx,Hy,range] = wgmode3d( 'geo_inj.in',0.014,1,1,2.5);
savetfsfinj( 'tfsfex.set', Ex, Hy, range, range );
savetfsfinj( 'tfsfey.set', Ey, Hx, range, range );
quit;
