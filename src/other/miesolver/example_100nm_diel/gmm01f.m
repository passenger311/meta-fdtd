clear all;

min_lambda = 400; %in nm
max_lambda = 800;
position = [0. 0. 0.];
%min_radius = 1;
%max_radius = 10;


frequ_factor = 2.99792458e5;  % change from frequency in THz (c|=1) to inverse wavelength in 1/nm (c=1)

% Drude-Lorentz material in THz
eps_inf = 5.9673;
omegaD = 2113.6;      % Drude plasma frequency [2 pi c]
gammaD = 15.92;  % Drude damping constant [1/dt]
omegaL = 650.07;      % Lorentzian plasma frequency [2 pi c]
gammaL = 104.86/2;   % Lorentzian damping constant [1/dt]
deltaeps = 1.09;

%for sqradius = min_radius : 1 : max_radius
  dirname = sprintf('%s','data');
  mkdir(dirname);

  for j = min_lambda:100:max_lambda

 %   In1 = (frequ_factor)/j;
 %   n = sqrt(eps_inf - omegaD^2./(In1*In1-i*In1*gammaD) - deltaeps*omegaL^2/(In1*In1-i*2.*gammaL*In1-omegaL^2));

    fid = fopen('input.in','w');
    fprintf(fid,'%e\n',j*1e-9);
    fprintf(fid,'%i\n',1);
    fprintf(fid,'%e,%e,%e,%e,%e,%e\n',position,100e-9,3.47,0.0);
    fclose(fid);
    !gmm01f
    !sed -i 1,+15d Mgmm01f.out
%    destination = sprintf('%s/%i.out',dirname,j);
%    movefile('Mgmm01f.out',destination);
    fid=fopen('Mgmm01f.out','r');
    fid_wr=fopen('Mgmm01f.dat','w');
    while feof(fid) == 0
       tline = fgetl(fid);
       tmp = sscanf(tline,'%f %e %e %e %e %e %e %e');
       fprintf(fid_wr,'%f %e %e\n', tmp(1),tmp(7)/pi, tmp(4)/pi);
    end
    fclose(fid); fclose(fid_wr);
					 


    destination = sprintf('%s/%i.out','data',j);
    movefile('Mgmm01f.dat',destination);

  end
%  cd(dirname);
%  addpath ..;
%  backscattering;
%  cd ..;
%end
