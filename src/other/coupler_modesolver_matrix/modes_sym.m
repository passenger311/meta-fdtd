clear all;

addpath tools


lambda = 1550;         % wavelength (nm)
neffguess = 2.29726;
nmodes = 1;            % number of modes to compute
dx = 10;               % grid size (x)
dy = dx;

% --- geometric quantities used to set up geometry

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
metaltop = 0;
metalside = 0;


nsubstrate = 1.444;%3.47;%    % SiO2 substrate
nsubsubbase = 3.47;%1.444
nsubbase = 1.444;%3.47;%1.444;%
nbase = 3.47;%1.444;%3.47;%          % Silicon base
ncore = 3.47;          % Silicon core
%ngap = 1.444;          % Gap
nsuperstrate = 1.00;
%metaltop=1;
%metalenclose=0;
%metalside=1;
epsr = -125;            % real part Metal on sides (on top)
epsi = 3;               % imaginary part Metal on sides (on top)

hsubstrate = 800;      % SiO2 substrate (nm)
%hbase = 0;             % Silicon base (nm)
%hcore = 220;           % Silicon core (nm)
%hgap = 0;           % Top gap (nm)
%hmetal = 50;
if hgap == 0; ngap = 1; end;
hsuperstrate = 800;    % Superstrate (nm)
wcore = 400;           % waveguide full-width (nm)
%wgap = 50;             % gap between waveguide and metal on sides (nm)
%wmetalenclose = 50;
%wmetal = 100;
wside = 600;            % metal on side of waveguide (nm)

% ---- setup geometry

fprintf (1,'generating index mesh...\n');
h = [hsubstrate,50,hsubsubbase,hsubbase,hbase,hcore,hgap,hmetal,hsuperstrate];
w = [wcore/2,wgap,wside];%wgap,wmetal,wside];
n = zeros(size(w,2),size(h,2));
n(:,1)=nsubstrate;
n(:,2)=nsubsubbase;
n(:,3)=nsubsubbase;
n(:,4)=nsubbase;
n(:,5)=nbase;
n(:,6)=nsuperstrate;
n(:,7)=nsuperstrate;
n(:,8)=nsuperstrate;
n(:,9)=nsuperstrate;
n(1,6)=ncore;
n(1:2,8)=sqrt(epsr+i*epsi);
n(2,6)=sqrt(epsr+i*epsi);

%n(3,4)=nbase;
n(3,2)=nsubstrate;

% ---- setting up the mesh

[x,y,xc,yc,nx,ny,eps] = waveguidemesh(n.*n,h,w,dx,dy);

% ---- Now we stretch out the mesh at the boundaries:

%[x,y,xc,yc,dx,dy] = stretchmesh(x,y,[round(hsuperstrate/dy/2),round(hsubstrate/dy/2),round(wside/dx/2),0],[4,4,4,0]);
[x,y,xc,yc,dx,dy] = stretchmesh(x,y,[round(hsuperstrate/dy/2),round(hsubstrate/dy/2),round(wside/dx/2),0],[1+j*2,1+j*2,4,0]);

% ---- Invoke mode calculator

j=1;
%for wcore = 200%:-40:0
  fprintf('\n\nCore thickness: %i\n',wcore);
  tmp=couplermode_sym(lambda,neffguess,nmodes,dx,dy,x,y,xc,yc,eps);
  neff(j,:)=tmp(:);
  fprintf(1,'neff = %7.5f -i*%7.5f\n',[real(neff(j,:));abs(imag(neff(j,:)))]);
%  j=j+1;
%end;
