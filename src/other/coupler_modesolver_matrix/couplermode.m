function [neff]=couplermode(lambda,wcore,wgap,wmetal,hsubsubbase,hsubbase,hbase,hcore,hgap,hmetal,ngap,neffguess,nmodes,metaltop,metalside,dx,dy)

addpath tools

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
%wcore = 400;           % waveguide full-width (nm)
%wgap = 50;             % gap between waveguide and metal on sides (nm)
%wmetalenclose = 50;
%wmetal = 100;
wside = 600;            % metal on side of waveguide (nm)

%dx = 4;                % grid size (x)
%dy = dx;               % grid size (y)

%lambda = 1550;         % wavelength (nm)
%nmodes = 1;            % number of modes to compute

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

[x,y,xc,yc,nx,ny,eps] = waveguidemesh(n,h,w,dx,dy);
% Now we stretch out the mesh at the boundaries:
%[x,y,xc,yc,dx,dy] = stretchmesh(x,y,[round(hsuperstrate/dy/2),round(hsubstrate/dy/2),round(wside/dx/2),0],[4,4,4,0]);
[x,y,xc,yc,dx,dy] = stretchmesh(x,y,[round(hsuperstrate/dy/2),round(hsubstrate/dy/2),round(wside/dx/2),0],[1+j*2,1+j*2,4,0]);

figure(2)
imagemode(xc,yc,log(real(eps))/log(max(max(abs(eps)))),'Epsilon structure');


fprintf (1,'solving for eigenmodes...'); t = cputime;
[Hx,Hy,neff] = wgmodes (lambda, ncore, nmodes, dx, dy, eps, '000S');
fprintf (1,'done (cputime = %7.3f)\n', cputime-t);



if nmodes==1,
   fprintf (1,'post-processing...'); t = cputime;
   [Hz,Ex,Ey,Ez] = postprocess (lambda, neff, Hx, Hy, dx, dy, eps, '000S');
   fprintf (1,'done (cputime = %7.3f)\n', cputime-t);

   ii = 1;
   colormap(jet(256));
   hn = max(abs([Hy(:);Hx(:)]));
   %en = max(abs([Ey(:);Ex(:)]));
   en = hn/neff;
   %en = hn/ncore;
   figure(1);
   subplot(231);
   imagemode(x,y,Hx/hn,sprintf('Hx (mode %d)',ii));
   hold on;
   contourmode(x,y,Hx/hn,1,3,45,sprintf('Hx (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(234);
   imagemode(x,y,Hy/hn,sprintf('Hy (mode %d)',ii));
   hold on;
   contourmode(x,y,Hy/hn,1,3,45,sprintf('Hy (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(233);
   imagemode(x,y,Hz/hn,sprintf('Hz (mode %d)',ii));
   hold on;
   contourmode(x,y,Hz/hn,1,3,45,sprintf('Hz (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(235);
   imagemode(xc,yc,Ex/en,sprintf('Ex (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ex/en,1,3,60,sprintf('Ex (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(232);
   imagemode(xc,yc,Ey/en,sprintf('Ey (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ey/en,1,3,60,sprintf('Ey (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(236);
   imagemode(xc,yc,Ez/en,sprintf('Ez (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ez/en,1,3,60,sprintf('Ez (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);
   figure(2)
   imagemode(xc,yc,log(real(eps))/log(max(max(abs(eps)))),'Epsilon structure');
end;
