% Parameters taken from Vial et al Phys. Rev. B, 2005, 71, 085416

clear all;

filename = sprintf('Au_eps_real.dat');
fid = fopen(filename,'w');
filename = sprintf('Au_eps_imag.dat');
fid2 = fopen(filename,'w');

eps_inf = 5.9673;
o_D = 2113.6e12*2*pi;
g_D = 15.92e12*2*pi;
de_L = 1.09;
o_L = 650.07e12*2*pi;
g_L = 104.86e12/2*2*pi;

for j=100:2:1550
  omeg = 299792458/(j*1e-9)*(2*pi);

  eps = eps_inf - o_D^2/(omeg^2-i*g_D*omeg) - de_L*o_L^2/(omeg^2-i*2*g_L*omeg-o_L^2);
  eps0 = eps_inf;
  eps1 = - o_D^2/(omeg^2-i*g_D*omeg);
  eps2 = - de_L*o_L^2/(omeg^2-i*2*g_L*omeg-o_L^2);
  fprintf(fid,'%e %e %e %e %e\n', j, real(eps), real(eps0), real(eps1), real(eps2));
  fprintf(fid2,'%e %e %e %e %e\n', j, imag(eps), imag(eps0), imag(eps1), imag(eps2));
end
fclose(fid); fclose(fid2);
exit
