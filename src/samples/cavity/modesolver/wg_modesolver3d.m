fid = fopen('lua.matlab', 'r');
tline = fgetl(fid); neff_guess = sscanf(tline,'%e');
tline = fgetl(fid); lambda_inv = sscanf(tline,'%e');
tline = fgetl(fid); blob_out = sscanf(tline,'%i %i %i %i %i %i %i %i %i');
tline = fgetl(fid); fileName_in = sscanf(tline,'%s');
tline = fgetl(fid); fileName_out = sscanf(tline,'%s');
fclose(fid);

lambda=1/lambda_inv;
nmodes = 1;         % number of modes to compute

addpath /home/atiraid/atitacpg/sw00100/Documents/MATLAB/modified_modesolver/tools

%start solver
[Hx,Hy,Hz,Ex,Ey,Ez,neff,blob] = wgmodes (fileName_in, lambda, neff_guess, nmodes, '0000');

fprintf(1,'neff = %15.12f\n',neff);


%save fields
fprintf('\nsaving fields in %s...\n\n', fileName_out)
saveFields (fileName_out,Ex,Hy,blob,blob_out);


%plot fields
ii = 1;
colormap(jet(256));
hn = max(abs(Hy(:)));
en = hn/(neff_guess);

dx = blob(3);
nrow = ceil((blob(2)-blob(1))/dx+.5);
dy = blob(6);
ncolumn = ceil((blob(5)-blob(4))/dy+.5);
xc = zeros(nrow,1); xc(:,1) = blob(1):dx:blob(2);
yc = zeros(1,ncolumn); yc(1,:) = blob(4):dy:blob(5);


%open file with eps values again to print 
fid = fopen(fileName_in, 'r');
if (fid==-1)
  fprintf('\nerror: cannot open %s for plot.\n', fileName_in);
end
fprintf(1,'reading file %s for plots...\n', fileName_in);

while feof(fid) == 0
   tline = fgetl(fid);
   if (length(tline)==0)
     continue;
   end
   if (tline(1)=='(')|(tline(2)=='!')
     continue;
   end
   [blob] = sscanf(tline,'%i %i %i %i %i %i %i %i %i');
   break;
end;

dx = blob(3);
nrow = ceil((blob(2)-blob(1))/dx+.5);
dy = blob(6);
ncolumn = ceil((blob(5)-blob(4))/dy+.5);
%read in file
epsxx=zeros(nrow,ncolumn);epsyy=zeros(nrow,ncolumn);
epszz=zeros(nrow,ncolumn);
row=1;column=1;

while feof(fid) == 0
   tline = fgetl(fid);
   if (length(tline)==0)
     continue;
   end
   if (tline(1)==')')
     continue;
   end
   tmp = sscanf(tline,'%e %e %e');
   epsxx(row,column)=tmp(1);
   epsyy(row,column)=tmp(2);
   epszz(row,column)=tmp(3);
   row = row+1;
   if(row>nrow);
     column=column+1;
     row=1;
   end;
end
fclose(fid);


fprintf('plotting results...\n');

figure(1);
clf(1);
surfc(sqrt(epsxx));
view(0,90);
shading interp;
colorbar;
axis equal;
axis tight;

print(gcf, '-dpng', '-r300', 'waveguide-structure.png')

figure(3);
clf(3);
subplot(321);
imagemode(xc,yc,Hx/hn,sprintf('Hx (mode %d)',ii));
hold on;
contourmode(xc,yc,Hx/hn,1,3,45,sprintf('Hx (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

subplot(322);
imagemode(xc,yc,Hy/hn,sprintf('Hy (mode %d)',ii));
hold on;
contourmode(xc,yc,Hy/hn,1,3,45,sprintf('Hy (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

subplot(325);
imagemode(xc,yc,Hz/hn,sprintf('Hz (mode %d)',ii));
hold on;
contourmode(xc,yc,Hz/hn,1,3,45,sprintf('Hz (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);


%fprintf (1,'post-processing...'); t = cputime;
%[Hz,Ex,Ey,Ez] = postprocess (lambda, neff, Hx, Hy, dx, dy, eps, '0000');
%fprintf (1,'done (cputime = %7.3f)\n', cputime-t);


subplot(324);
imagemode(xc,yc,Ex/en,sprintf('Ex (mode %d)',ii));
hold on;
contourmode(xc,yc,Ex/en,1,3,60,sprintf('Ex (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

subplot(323);
imagemode(xc,yc,Ey/en,sprintf('Ey (mode %d)',ii));
hold on;
contourmode(xc,yc,Ey/en,1,3,60,sprintf('Ey (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

subplot(326);
imagemode(xc,yc,Ez/en,sprintf('Ez (mode %d)',ii));
hold on;
contourmode(xc,yc,Ez/en,1,3,60,sprintf('Ez (mode %d)',ii));
hold off;
v = xlim();
%line(v,[h1,h1]);
%line([0,w/2,w/2],[h1+h2,h1+h2,h1]);
print(gcf, '-dpng', '-r300', 'waveguide-fields.png')


fprintf('finished.\n');
exit
