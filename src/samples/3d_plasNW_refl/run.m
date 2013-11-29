[Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( 'geo_injection',0.0052631578947368,1,1,2.9,0);
wgsave3k( 'tfsfex_re.in', real(Ex), real(Hy), range, range );
wgsave3k( 'tfsfex_im.in', imag(Ex), imag(Hy), range, range );
wgsave3k( 'tfsfey_re.in', real(Ey), -real(Hx), range, range );
wgsave3k( 'tfsfey_im.in', imag(Ey), -imag(Hx), range, range );
wgsave3k3( sprintf('Emode%2.2f_re.set',neff), real(Ex), real(Ey), real(Ez), range, range );
wgsave3k3( sprintf('Emode%2.2f_im.set',neff), imag(Ex), imag(Ey), imag(Ez), range, range );
wgsave3k3( sprintf('Hmode%2.2f_re.set',neff), real(Hx), real(Hy), real(Hz), range, range );
wgsave3k3( sprintf('Hmode%2.2f_im.set',neff), imag(Hx), imag(Hy), imag(Hz), range, range );
quit;
