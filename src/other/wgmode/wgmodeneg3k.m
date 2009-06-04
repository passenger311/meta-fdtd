function [ Ex,Ey,Hx,Hy,range] = wgmodeneg3k( epsfile, lambdainv, hx, hy, neffguess )
%  Detailed explanation goes here

om = 2*pi*lambdainv;

betaguess = (neffguess * om)^2;

% --- load .set file -> eps(1,l) ; l is the linear index

fprintf('loading file ...\n');

fid = fopen(epsfile,'r');

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

%l = 0;
%while feof(fid) == 0
%   tline = fgetl(fid);
%   %if (length(tline)==0)
%   %  continue;
%   %end
%   if (tline(1)==')')
%     continue;
%   end
%   l = l + 1;
%   tmp = sscanf(tline,'%e %e %e');
%   leps(l,1) = tmp(1);
%   leps(l,2) = tmp(2);
%   leps(l,3) = tmp(3);
%end 
tmp = textscan(fid, '%n %n %n %n %n %n');
l = size(tmp{1},1);
leps = zeros(l,size(tmp,2)/2);
lmu = zeros(l,size(tmp,2)/2);
for j = 1 : size(tmp,2)/2
   leps(:,j) = tmp{j};
   lmu(:,j) = tmp{size(tmp,2)/2+j};
end

%DELETE THE FOLLOWING LINE!!!!!!
%lmu=leps;

if ( l ~= sl )
 error('error')
end
fclose(fid);

if ( si<1 || sj<1 )
 error('error')
end

leps = [leps(1:si,:);leps(:,:);leps(sl-si+1:sl,:)]; %add a row at the beginning and one at the end
%leps(sl+1:sl+si,1) = leps(sl,1); % add a row
%leps(sl+1:sl+si,2) = leps(sl,2); % add a row
%leps(sl+1:sl+si,3) = leps(sl,3); % add a row 
lmu = [lmu(1:si,:);lmu(:,:);lmu(sl-si+1:sl,:)]; %add a row at the beginning and one at the end
%lmu(sl+1:sl+si,1) = lmu(sl,1); % add a row
%lmu(sl+1:sl+si,2) = lmu(sl,2); % add a row
%lmu(sl+1:sl+si,3) = lmu(sl,3); % add a row

% --- plot averaged refractive index

fprintf('plotting epsilon ...\n');

eps = zeros(sj,si);
mu = zeros(sj,si);

nmax = 0;
for j=1:sj
    for i=1:si
        m=j*si+i;
        eps(j,i) = leps(m,3); %0.25*(leps(l,1)+leps(l-1,1)+leps(l,2)+leps(l-si,2));
	mu(j,i) = lmu(m,3);
        if real(sqrt(eps(j,i)*mu(j,i))) > nmax 
           nmax = real(sqrt(eps(j,i)*mu(j,i)));
        end
    end
end

fprintf(' -> nmax = %f\n', nmax);

%eps(1,:) = eps(2,:);
%eps(:,1) = eps(:,2);
%mu(1,:) = mu(2,:);
%mu(:,1) = mu(:,2);

subplot(3,2,1)
hold on;
title('Refractive index profile n(x,y)')
xlabel('x [a.u.]')
ylabel('y [a.u.]')
surfc(real(sqrt(eps.*mu)));
axis tight;
shading interp;
colorbar;
view(0,90);
hold off;

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
    for i=1:si
        l=(j-1)*si+i; m = l + si;  
        % --- block 1 (Hx-Hx)
        A(1,l) = leps(m,2)/leps(m,3)/(hy*hy);
        C(1,l) = lmu(m-1,1)/lmu(m-1,3)/(hx*hx);
        D(1,l) = -lmu(m,1)/(hx*hx)*(1./lmu(m,3)+1./lmu(m-1,3)) + lmu(m,1)*leps(m,2)*om*om - leps(m,2)/(hy*hy)*(1./leps(m+si,3)+1./leps(m,3));
        E(1,l) = lmu(m+1,1)/lmu(m,3)/(hx*hx);
        G(1,l) = leps(m,2)/leps(m+si,3)/(hy*hy); 
        % --- block 2 (Hx-Hy)
        C(2,l) = 1./(hx*hy)*(lmu(m-1,2)/lmu(m-1,3)-leps(m,2)/leps(m,3));
        D(2,l) = -1./(hx*hy)*(lmu(m,2)/lmu(m,3)-leps(m,2)/leps(m,3));
        F(2,l) = -1./(hx*hy)*(lmu(m+si-1,2)/lmu(m-1,3)-leps(m,2)/leps(m+si,3));
        G(2,l) = 1./(hx*hy)*(lmu(m+si,2)/lmu(m,3)-leps(m,2)/leps(m+si,3));
        % --- block 3 (Hy-Hy)
        A(3,l) = lmu(m-si,2)/lmu(m-si,3)/(hy*hy);
        C(3,l) = leps(m,1)/leps(m,3)/(hx*hx);
        D(3,l) = -lmu(m,2)/(hy*hy)*(1./lmu(m,3)+1./lmu(m-si,3)) + lmu(m,2)*leps(m,1)*om*om - leps(m,1)/(hx*hx)*(1./leps(m+1,3)+1./leps(m,3));
        E(3,l) = leps(m,1)/leps(m+1,3)/(hx*hx);
        G(3,l) = lmu(m+si,2)/lmu(m,3)/(hy*hy);
        % --- block 4 (Hy-Hx)
        A(4,l) = 1./(hx*hy)*(lmu(m-si,1)/lmu(m-si,3)-leps(m,1)/leps(m,3));
        B(4,l) = -1./(hx*hy)*(lmu(m-si+1,1)/lmu(m-si,3)-leps(m,1)/leps(m+1,3));
        D(4,l) = -1./(hx*hy)*(lmu(m,1)/lmu(m,3)-leps(m,1)/leps(m,3));
        E(4,l) = 1./(hx*hy)*(lmu(m+1,1)/lmu(m,3)-leps(m,1)/leps(m+1,3));       
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

Hx = zeros(sj,si);
Hy = zeros(sj,si);

for j=1:sj
    for i=1:si
        l=(j-1)*si+i;
        Hx(j,i) = EVec(l,1);
        Hy(j,i) = EVec(l+sl,1);
    end
end

% --- calculate remaining field components

fprintf('calculating E field components ...\n');

Hz = zeros(sj,si);
Ex = zeros(sj,si);
Ey = zeros(sj,si);

for j=1:sj-1
    for i=1:si-1
        m=j*si+i;
        Hz(j,i) = -1/keff/lmu(m,3)*( (lmu(m+1,1)*Hx(j,i+1)-lmu(m,1)*Hx(j,i))/hx + (lmu(m+si,2)*Hy(j+1,i)-lmu(m,2)*Hy(j,i))/hy); % = Hz / i
    end
end

Hz(sj,:) = Hz(sj-1,:);
Hz(:,si) = Hz(:,si-1);

emax = 0;
for j=2:sj
    for i=2:si
        m=j*si+i;
        Ex(j,i) = 1./(leps(m,1)*om) * ( (Hz(j,i)-Hz(j-1,i))/hy + keff*Hy(j,i) ); 
        Ey(j,i) = - 1./(leps(m,2)*om) * ( (Hz(j,i)-Hz(j,i-1))/hx + keff*Hx(j,i) ); 
        if abs(Ex(j,i)) > emax 
            emax = abs(Ex(j,i));
        end
        if abs(Ey(j,i)) > emax 
            emax = abs(Ey(j,i));
        end
    end
end

Ex(1,:) = Ex(2,:);
Ex(:,1) = Ex(:,2);
Ey(1,:) = Ey(2,:);
Ey(:,1) = Ey(:,2);

Ex = real(Ex)/emax;
Ey = real(Ey)/emax;
Hx = real(Hx)/emax;
Hy = real(Hy)/emax;

% --- plot fields

fprintf('plotting fields ...\n');

colormap(jet(256));
cmax = size(colormap,1)-1;
v = (0:-3:-45)';

subplot(3,2,3);

hold on;
%contour(Ex',25);
title(sprintf('Ex component @ neff = %f',neff));
ylabel('x [a.u.]')
xlabel('y [a.u.]')
modebmp = uint8(abs(cmax*Ex));
image(modebmp);
set(gca,'YDir','normal');
axis equal;
axis tight;
contour(20*log10(abs(Ex)),v);
hold off;

subplot(3,2,4);

hold on;
%contour(Ey',25);
title(sprintf('Ey component @ neff = %f',neff));
ylabel('x [a.u.]')
xlabel('y [a.u.]')
modebmp = uint8(abs(cmax*Ey));
image(modebmp);
set(gca,'YDir','normal');
axis equal;
axis tight;
contour(20*log10(abs(Ey)),v);
hold off;

subplot(3,2,6);

hold on;
%contour(Hx',25);
title(sprintf('Hx component @ neff = %f',neff));
ylabel('x [a.u.]')
xlabel('y [a.u.]')
modebmp = uint8(abs(cmax*Hx));
image(modebmp);
set(gca,'YDir','normal');
axis equal;
axis tight;
contour(20*log10(abs(Hx)),v);
hold off;

subplot(3,2,5);
hold on;
title(sprintf('Hy component @ neff = %f',neff));
ylabel('x [a.u.]')
xlabel('y [a.u.]')
modebmp = uint8(abs(cmax*Hy));
image(modebmp);
set(gca,'YDir','normal');
axis equal;
axis tight;
contour(20*log10(abs(Hy)),v);
hold off;

fprintf('done.\n');






