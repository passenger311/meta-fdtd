[ Ex,Ey,Hx,Hy,range] = wgmode3d( 'geo_inj.in',0.0096061479346782,1,1,2.5);
savetfsfinj( 'tfsfex.set', Ex, Hy, range, range );
savetfsfinj( 'tfsfey.set', Ey, Hx, range, range );
quit;
