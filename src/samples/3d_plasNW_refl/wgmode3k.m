function [ Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( epsfile, lambdainv, hx, hy, neffguess, fig_on )
%  Detailed explanation goes here

om = 2*pi*lambdainv;

betaguess = (neffguess * om)^2;

% --- load .set file -> eps(1,l) ; l is the linear index

fprintf('loading file ...\n');

fid = fopen(sprintf('%s_re.in',epsfile),'r');

while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' ) 
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end


ri = [ range(1) range(2) ]; % range i
rj = [ range(4) range(5) ]; % range j

si = ri(2)-ri(1)+1; % length i
sj = rj(2)-rj(1)+1; % length j

sl = si*sj;

leps = zeros(sl,3);

l = 0;
while feof(fid) == 0
   tline = fgetl(fid);
   %if (length(tline)==0)
   %  continue;
   %end
   if (tline(1)==')')
     continue;
   end
   l = l + 1;
   tmp = sscanf(tline,'%e %e %e');
   leps(l,1) = tmp(1);
   leps(l,2) = tmp(2);
   leps(l,3) = tmp(3);
end 
if ( l ~= sl )
 error('error')
end
fclose(fid);

if ( si<1 || sj<1 )
 error('error')
end


fid = fopen(sprintf('%s_im.in',epsfile),'r');

while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' ) 
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end


l = 0;
while feof(fid) == 0
   tline = fgetl(fid);
   %if (length(tline)==0)
   %  continue;
   %end
   if (tline(1)==')')
     continue;
   end
   l = l + 1;
   tmp = sscanf(tline,'%e %e %e');
   leps(l,1) = leps(l,1)+i*tmp(1);
   leps(l,2) = leps(l,2)+i*tmp(2);
   leps(l,3) = leps(l,3)+i*tmp(3);
end 
if ( l ~= sl )
 error('error')
end
fclose(fid);


leps(sl+1:sl+si,1) = leps(sl,1); % add a row
leps(sl+1:sl+si,2) = leps(sl,2); % add a row
leps(sl+1:sl+si,3) = leps(sl,3); % add a row 

fprintf('plotting epsilon ...\n');

eps = zeros(sj,si);

nmax = 0;
for j=2:sj
    for ii=2:si
        l=(j-1)*si+ii;
        eps(j,ii) = leps(l,3); %0.25*(leps(l,1)+leps(l-1,1)+leps(l,2)+leps(l-si,2));
        if abs(sqrt(eps(j,ii))) > nmax 
           nmax = abs(sqrt(eps(j,ii)));
        end
    end
end

fprintf(' -> nmax = %f\n', nmax);

eps(1,:) = eps(2,:);
eps(:,1) = eps(:,2);

if (fig_on==1);
clf
figure(1)
scrsz = get(0,'ScreenSize');
%figure('Position',[100 1 scrsz(3)/3 scrsz(4)-249])
%set(gcf,'Position',[100 1 scrsz(3)-200 scrsz(4)-249])
%set(gcf,'Position',[0 0 scrsz(3) scrsz(4)])
set(gcf,'Position',[0 0 1000 600])

subplot(3,3,1)
hold on;
title('Refractive index profile n(x,z)')
xlabel('x [a.u.]')
ylabel('z [a.u.]')
pcolor(real(eps));
axis tight;
shading interp;
colorbar;
%view(0,9);
hold off;
end

% --- calculate diagonals A,B,C,D,E,F for 4 BLOCKS

fprintf('calculating matrix elements ...\n');

% note: linear index l relates to i,j as follows:
% l + 1  -> i+1,j
% l - 1  -> i-1,j
% l + si -> i,j+1
% l - si -> i,j-1
% l + sl -> Hy

A = zeros(4,sl); % A -> i,j-1
B = zeros(4,sl); % B -> i+1,j-1
C = zeros(4,sl); % C -> i-1,j
D = zeros(4,sl); % D -> i,j
E = zeros(4,sl); % E -> i+1,j
F = zeros(4,sl); % F -> i-1,j+1
G = zeros(4,sl); % G -> i,j+1

% in real (i,j) coordinates
%  j
%  ^
%  +  F G
%  +  C D E
%  +    A B 
%  +--+-+-+---> i 

% matrix block arrangement:
%  1 2 
%  4 3

for j=1:sj
    for ii=1:si
        l=(j-1)*si+ii;   
        % --- block 1 (Hx-Hx)
        A(1,l) = leps(l,2)/leps(l,3)/(hy*hy);
        C(1,l) = 1./(hx*hx);
        D(1,l) = -2./(hx*hx) + leps(l,2)*om*om - leps(l,2)/(hy*hy)*(1./leps(l+si,3)+1./leps(l,3));
        E(1,l) = 1./(hx*hx);
        G(1,l) = leps(l,2)/leps(l+si,3)/(hy*hy); 
        % --- block 2 (Hx-Hy)
        C(2,l) = 1./(hx*hy)*(1-leps(l,2)/leps(l,3));
        D(2,l) = -C(2,l);
        F(2,l) = -1./(hx*hy)*(1-leps(l,2)/leps(l+si,3));
        G(2,l) = -F(2,l);
        % --- block 3 (Hy-Hy)
        A(3,l) = 1./(hy*hy);
        C(3,l) = leps(l,1)/leps(l,3)/(hx*hx);
        D(3,l) = -2./(hy*hy) + leps(l,1)*om*om - leps(l,1)/(hx*hx)*(1./leps(l+1,3)+1./leps(l,3));
        E(3,l) = leps(l,1)/leps(l+1,3)/(hx*hx);
        G(3,l) = 1./(hy*hy);
        % --- block 4 (Hy-Hx)
        A(4,l) = 1./(hx*hy)*(1-leps(l,1)/leps(l,3));
        B(4,l) = -1./(hx*hy)*(1-leps(l,1)/leps(l+1,3));
        D(4,l) = -A(4,l);
        E(4,l) = -B(4,l);       
    end
end


% --- assume zero derivative boundary condition

fprintf('fixing boundary elements ...\n');

% south edge
%D(:,1:si) = D(:,1:si)+ A(:,1:si);
%E(:,1:si) = E(:,1:si)+ B(:,1:si);
A(:,1:si) = 0.;
B(:,1:si) = 0.;

% north edge
o = (sj-1)*si; 
%D(:,o+1:o+si) = D(:,o+1:o+si)+ G(:,o+1:o+si);
%C(:,o+1:o+si) = C(:,o+1:o+si)+ F(:,o+1:o+si);
G(:,o+1:o+si) = 0.;
F(:,o+1:o+si) = 0.;

% west edge
%D(:,1:si:sl) = D(:,1:si:sl)+ C(:,1:si:sl);
%G(:,1:si:sl) = G(:,1:si:sl)+ F(:,1:si:sl);
C(:,1:si:sl) = 0.;
F(:,1:si:sl) = 0.;

% east edge
%D(:,si:si:sl) = D(:,si:si:sl)+ E(:,si:si:sl);
%A(:,si:si:sl) = A(:,si:si:sl)+ B(:,si:si:sl);
E(:,si:si:sl) = 0.;
B(:,si:si:sl) = 0.;

% --- define column indices for diagonals

Am = [1:sl] - si;
Bm = [1:sl] - si + 1;
Cm = [1:sl] - 1;
Dm = [1:sl];
Em = [1:sl] + 1;
Fm = [1:sl] + si - 1;
Gm = [1:sl] + si;

% --- build sparse matrix with indices (l,m) for (row,column)

fprintf('setting sparse matrix elements ...\n');

loff = [ 0, 0, sl, sl ]; % block row offset 
moff = [ 0, sl, sl, 0 ]; % block column offset

Lt = [];
Mt = [];
Vt = [];

l = [1:sl];
for b=1:4 
    lo = loff(b);
    mo = moff(b);
    Lt = [ Lt, l+lo, l+lo, l+lo, l+lo, l+lo, l+lo, l+lo];
    Mt = [ Mt, Am+mo, Bm+mo, Cm+mo, Dm+mo, Em+mo, Fm+mo, Gm+mo ];
    Vt = [ Vt, A(b,:),B(b,:),C(b,:),D(b,:),E(b,:),F(b,:),G(b,:) ];
end

c = 0;
for k = 1:length(Vt) % cleanup matrix
    v = Vt(k);
    if v ~= 0 
        c = c+1;
    end
end

L = zeros(1,c);
M = zeros(1,c);
V = zeros(1,c);

fprintf('wiping out zero values ...\n');

c = 0;
for k = 1:length(Vt) % cleanup matrix
    v = Vt(k);
    if v ~= 0
        c = c + 1;
        L(1,c) = Lt(k);
        M(1,c) = Mt(k);
        V(1,c) = v;
    end
end

fprintf('constructing sparse matrix ...\n');

S = sparse(L,M,V,sl*2,sl*2,length(V));

% --- run eigensolver

fprintf('running sparse solver ...\n');

opts.disp = 0;
[EVec,EVal] = eigs(S,[],1, betaguess,opts);

keff = sqrt(EVal(1));

neff = keff/om;
fprintf(' -> keff = %f\n', keff);
fprintf(' -> neff = %f\n', neff);

f_out=fopen('neff.lua','w');
fprintf(f_out,'nrefr = %f\n',neff);
fprintf(f_out,'nrefr_im = %f\n',imag(neff));
fclose(f_out);

Hx = zeros(sj,si);
Hy = zeros(sj,si);

for j=1:sj
    for ii=1:si
        l=(j-1)*si+ii;
        Hx(j,ii) = EVec(l,1);
        Hy(j,ii) = EVec(l+sl,1);
    end
end

% --- calculate remaining field components

fprintf('calculating E field components ...\n');

Hz = zeros(sj,si);
Ex = zeros(sj,si);
Ey = zeros(sj,si);

for j=1:sj-1
    for ii=1:si-1
        Hz(j,ii) = -i/keff*( (Hx(j,ii+1)-Hx(j,ii))/hx + (Hy(j+1,ii)-Hy(j,ii))/hy); % = Hz / ii
    end
end

Hz(sj,:) = Hz(sj-1,:);
Hz(:,si) = Hz(:,si-1);

emax = 0;
for j=2:sj
    for ii=2:si
        l=(j-1)*si+ii;
        Ex(j,ii) = 1./(leps(l,1)*om) * ( -i*(Hz(j,ii)-Hz(j-1,ii))/hy + keff*Hy(j,ii) );
        Ey(j,ii) = - 1./(leps(l,2)*om) * ( -i*(Hz(j,ii)-Hz(j,ii-1))/hx + keff*Hx(j,ii) );
        Ez(j,ii) = i./(leps(l,3)*om) * ( (Hx(j,ii)-Hx(j-1,ii))/hy - (Hy(j,ii)-Hy(j,ii-1))/hx );
        if abs(Ex(j,ii)) > abs(emax)
            emax = Ex(j,ii);
        end
        if abs(Ey(j,ii)) > abs(emax)
            emax = Ey(j,ii);
        end
    end
end

Ex(1,:) = Ex(2,:);
Ex(:,1) = Ex(:,2);
Ey(1,:) = Ey(2,:);
Ey(:,1) = Ey(:,2);
Ez(1,:) = Ez(2,:);
Ez(:,1) = Ez(:,2);

Ex = Ex/emax;
Ey = Ey/emax;
Ez = Ez/emax;
Hx = Hx/emax;
Hy = Hy/emax;
Hz = Hz/emax;

%for j=2:sj
%    for ii=2:si
%        l=(j-1)*si+ii;
%        Ex(j,ii) = 1./(leps(l,1)*om) * ( (Hz(j,ii)-Hz(j-1,ii))/hy + keff*Hy(j,ii) ); 
%        Ey(j,ii) = - 1./(leps(l,2)*om) * ( (Hz(j,ii)-Hz(j,ii-1))/hx + keff*Hx(j,ii) );
%        Ez(j,ii) = 1./(leps(l,3)*om) * ( (Hx(j,ii)-Hx(j-1,ii))/hy - (Hy(j,ii)-Hy(j,ii-1))/hx );
%    end
%end
%
%Ex(1,:) = Ex(2,:);
%Ex(:,1) = Ex(:,2);
%Ey(1,:) = Ey(2,:);
%Ey(:,1) = Ey(:,2);
%Ez(1,:) = Ez(2,:);
%Ez(:,1) = Ez(:,2);
%
%emax = sqrt(max(max((Ex.*Ex)'+(Ey.*Ey)'+(Ez.*Ez)')));%max(max(abs([Ex Ey Ez])))
%hmax = sqrt(max(max((Hx.*Hx)'+(Hy.*Hy)'+(Hz.*Hz)')));%max(max(abs([Hx Hy Hz])));
%
%Ex = Ex/emax;
%Ey = Ey/emax;
%Ez = Ez/emax;
%Hx = Hx/hmax;
%Hy = Hy/hmax;
%Hz = Hz/hmax;

% --- plot fields

if (fig_on==1);
fprintf('plotting fields ...\n');

subplot(3,3,2);
hold on;
%contour(Ex',25);
title(sprintf('El. field intensity @ neff = %5.3f',neff));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor((Ex.*conj(Ex))'+(Ey.*conj(Ey))'+(Ez.*conj(Ez))');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,4);
hold on;
%contour(Ex',25);
title(sprintf('Ex component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(real(Ex)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,5);
hold on;
%contour(Ey',25);
title(sprintf('Ey component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(real(Ey)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,6);
hold on;
%contour(Ey',25);
title(sprintf('Ez component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(real(Ez)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,3);
hold on;
%contour(Hx',25);
title(sprintf('Mag. field intensity @ neff = %5.3f',neff));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor((Hx.*conj(Hx))'+(Hy.*conj(Hy))'+(Hz.*conj(Hz))');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,7);
hold on;
%contour(Hx',25);
title(sprintf('Hx component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(real(Hx)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,8);
hold on;
title(sprintf('Hy component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
%contour(Hy',25);
pcolor(real(Hy)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,9);
hold on;
title(sprintf('Hz component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
%contour(Hy',25);
pcolor(real(Hz)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

colormap jet(512);

filename=sprintf('mode_%5.3f_i%5.3f_re.png',real(neff),imag(neff));
%filename2=sprintf('mode_%5.3f.fig',neff);
%saveas(figure(1), filename2);
print(gcf, '-dpng', '-r300', filename);

figure(2)
scrsz = get(0,'ScreenSize');
%figure('Position',[100 1 scrsz(3)/3 scrsz(4)-249])
%set(gcf,'Position',[100 1 scrsz(3)-200 scrsz(4)-249])
%set(gcf,'Position',[0 0 scrsz(3) scrsz(4)])
set(gcf,'Position',[0 0 1000 600])
subplot(3,3,1)
hold on;
title('Refractive index profile n(x,z)')
xlabel('x [a.u.]')
ylabel('z [a.u.]')
pcolor(real(eps));
axis tight;
shading interp;
colorbar;
%view(0,9);
hold off;

subplot(3,3,2);
hold on;
%contour(Ex',25);
title(sprintf('El. field intensity @ neff = %5.3f',neff));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor((Ex.*conj(Ex))'+(Ey.*conj(Ey))'+(Ez.*conj(Ez))');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,4);
hold on;
%contour(Ex',25);
title(sprintf('Ex component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(imag(Ex)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,5);
hold on;
%contour(Ey',25);
title(sprintf('Ey component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(imag(Ey)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,6);
hold on;
%contour(Ey',25);
title(sprintf('Ez component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(imag(Ez)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,3);
hold on;
%contour(Hx',25);
title(sprintf('Mag. field intensity @ neff = %5.3f',neff));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor((Hx.*conj(Hx))'+(Hy.*conj(Hy))'+(Hz.*conj(Hz))');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,7);
hold on;
%contour(Hx',25);
title(sprintf('Hx component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
pcolor(imag(Hx)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,8);
hold on;
title(sprintf('Hy component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
%contour(Hy',25);
pcolor(imag(Hy)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

subplot(3,3,9);
hold on;
title(sprintf('Hz component'));
ylabel('x [a.u.]')
xlabel('z [a.u.]')
%contour(Hy',25);
pcolor(imag(Hz)');
axis tight;
shading interp;
colorbar;
view(90,-90);
hold off;

colormap jet(512);

filename=sprintf('mode_%5.3f_i%5.3f_im.png',real(neff),imag(neff));
%filename2=sprintf('mode_%5.3f.fig',neff);
%saveas(figure(1), filename2);
print(gcf, '-dpng', '-r300', filename);

close(figure);
end

fprintf('done.\n');






