% Parameters taken from Supplementary Information of McMahon et al, J. Phys. Chem. C, 2009, 113, 2731-2735
% Only publication contained wrong parameters (June 2009), this could have been changed by now
% The materials are Carbon, Silver measured by Johnson and Christy and Silver measured by Lynch and Hunter.

clear all;
for k=1:3

  f = {'C','Ag_JC','Ag_LH'};
  filename = sprintf('%s_eps_real.dat',f{k});
  fid = fopen(filename,'w');
  filename = sprintf('%s_eps_imag.dat',f{k})
  fid2 = fopen(filename,'w');
  if  k == 1%'carbon'
    eps_inf = 2.554;
    o_D = 2.14977e15;
    g_D = 9.47568e15;
    de_L1 = 3.56779;
    o_L1 = 6.82e15;
    g_L1 = 4.79937e15;
    de_L2 = 2.51021;
    o_L2 = 2.83344e15;
    g_L2 = 1.82616e15;
  elseif k == 2 %'silver'
    eps_inf = 1.17152;
    o_D = 1.39604e16;
    g_D = 12.6126;
    de_L1 = 2.23994;
    o_L1 = 8.25718e15;
    g_L1 = 1.95614e14;
    de_L2 = 0.222651;
    o_L2 = 3.05707e15;
    g_L2 = 8.52675e14;
  elseif k == 3%'silver2'
    eps_inf = 2.3646;
    o_D = 1.32749e16;
    g_D = 1.13778e14;
    de_L1 = 0.31506;
    o_L1 = 6.6547e15;
    g_L1 = 4.25395e14;
    de_L2= 0.86804;
    o_L2 = 7.87437e15;
    g_L2 = 8.32863e14;
  end


  for j=100:2:1550
    omeg = 299792458/(j*1e-9)*(2*pi);

    eps = eps_inf - o_D^2/(omeg^2-i*g_D*omeg) - de_L1*o_L1^2/(omeg^2-i*2*g_L1*omeg-o_L1^2) - de_L2*o_L2^2/(omeg^2-i*2*g_L2*omeg-o_L2^2);
    eps0 = eps_inf;
    eps1 = - o_D^2/(omeg^2-i*g_D*omeg);
    eps2 = - de_L1*o_L1^2/(omeg^2-i*2*g_L1*omeg-o_L1^2);
    eps3 = - de_L2*o_L2^2/(omeg^2-i*2*g_L2*omeg-o_L2^2);
    fprintf(fid,'%e %e %e %e %e %e\n', j, real(eps), real(eps0), real(eps1), real(eps2), real(eps3));
    fprintf(fid2,'%e %e %e %e %e %e\n', j, imag(eps), imag(eps0), imag(eps1), imag(eps2), imag(eps3));
  end
  fclose(fid); fclose(fid2);
end;
exit
