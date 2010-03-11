clear all;
lambda = 1550;         % wavelength (nm)
%wcore = 420;           % core width (nm)
wgap = 50;             % gap between waveguide and metal on sides (nm)
wmetal = 600;           % metal width on side of gap/waveguide (nm)
hsubsubbase=50;
hsubbase = 50;
hbase = 50;             % Silicon base (nm)
hcore = 200;           % Silicon core (nm)
hgap = 0;             % Top gap (nm)
hmetal = 50;           % Metal cover (nm)
ngap = 1.444;          % Gap
neffguess = 2.29726;
nmodes = 1;            % number of modes to compute
metaltop = 0;
metalside = 0;
dx = 10;               % grid size (x)
dy = dx;

j=1;
for wcore = 200%:-40:0
  fprintf('\n\nCore thickness: %i\n',wcore);
  tmp=couplermode(1550,wcore,wgap,wmetal,hsubsubbase,hsubbase,hbase,hcore,hgap,hmetal,ngap,neffguess,nmodes,metaltop,metalside,dx,dy);
  neff(j,:)=tmp(:);
  fprintf(1,'neff = %7.5f -i*%7.5f\n',[real(neff(j,:));abs(imag(neff(j,:)))]);
  j=j+1;
end;
