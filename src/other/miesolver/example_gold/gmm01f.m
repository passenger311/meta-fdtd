clear all;

min_lambda = 500; %in nm
max_lambda = 1000;
position = [0. 0. 0.];
min_radius = 1;
max_radius = 10;


frequ_factor = 2.99792458e5;  % change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

% Drude-Lorentz material in THz
eps_inf = 5.9673;
omegaD = 2113.6;      % Drude plasma frequency [2 pi c]
gammaD = 15.92;  % Drude damping constant [1/dt]
omegaL = 650.07;      % Lorentzian plasma frequency [2 pi c]
gammaL = 104.86/2;   % Lorentzian damping constant [1/dt]
deltaeps = 1.09;

for sqradius = min_radius : 1 : max_radius
  dirname = sprintf('data%ie-9',sqradius^2*4);
  mkdir(dirname);

  for j = min_lambda:max_lambda

    In1 = (frequ_factor)/j;
    n = sqrt(eps_inf - omegaD^2./(In1*In1-i*In1*gammaD) - deltaeps*omegaL^2/(In1*In1-i*2.*gammaL*In1-omegaL^2));

    fid = fopen('input.in','w');
    fprintf(fid,'%e\n',j*1e-9);
    fprintf(fid,'%i\n',1);
    fprintf(fid,'%e,%e,%e,%e,%e,%e\n',position,sqradius^2*4e-9,real(n),imag(n));
    fclose(fid);
    !gmm01f
    !sed -i 1,+15d Mgmm01f.out
    destination = sprintf('%s/%i.out',dirname,j);
    movefile('Mgmm01f.out',destination);

  end
  cd(dirname);
  addpath ..;
  backscattering;
  cd ..;
end
