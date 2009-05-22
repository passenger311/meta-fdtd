function [ Ex,Ey,Hx,Hy,range] = wgmode2j( epsfile, lambdainv, hx, hy, neffguess )
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

sl = si;

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

if ( si<1 )
 error('error')
end

leps(sl+1:sl+1,1) = leps(sl,1); % add a value
leps(sl+1:sl+1,2) = leps(sl,2); % add a value
leps(sl+1:sl+1,3) = leps(sl,3); % add a value

% --- plot averaged refractive index

fprintf('plotting epsilon ...\n');

eps = zeros(1,si);

nmax = 0;
for i=2:si
   l=i;
   eps(i) = leps(l,3); %0.25*(leps(l,1)+leps(l-1,1)+leps(l,2)+leps(l-si,2));
   if sqrt(eps(i)) > nmax 
     nmax = sqrt(eps(i));
   end
end

fprintf(' -> nmax = %f\n', nmax);

eps(1) = eps(2);

subplot(3,2,1)
hold on;
title('Refractive index profile n(x)')
xlabel('x [a.u.]')
plot(sqrt(eps));
axis tight;
hold off;

% --- calculate diagonals A,B,C,D,E,F for 4 BLOCKS

fprintf('calculating matrix elements ...\n');

% note: linear index l relates to i,j as follows:
% l + 1  -> i+1,j
% l - 1  -> i-1,j
% l + si -> i,j+1
% l - si -> i,j-1
% l + sl -> Hy

C = zeros(4,sl); % C -> i-1,j
D = zeros(4,sl); % D -> i,j
E = zeros(4,sl); % E -> i+1,j

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


    for i=1:si
        l=i;   
        % --- block 1 (Hx-Hx)
        C(1,l) = 1./(hx*hx);
        D(1,l) = -2./(hx*hx) + leps(l,2)*om*om;
        E(1,l) = 1./(hx*hx);
        % --- block 3 (Hy-Hy)
        C(3,l) = leps(l,1)/leps(l,3)/(hx*hx);
        D(3,l) = leps(l,1)*om*om - leps(l,1)/(hx*hx)*(1./leps(l+1,3)+1./leps(l,3));
        E(3,l) = leps(l,1)/leps(l+1,3)/(hx*hx);    
    end



% --- assume zero derivative boundary condition

fprintf('fixing boundary elements ...\n');

% west edge
%D(:,1:si:sl) = D(:,1:si:sl)+ C(:,1:si:sl);
%G(:,1:si:sl) = G(:,1:si:sl)+ F(:,1:si:sl);
C(:,1) = 0.;

% east edge
%D(:,si:si:sl) = D(:,si:si:sl)+ E(:,si:si:sl);
%A(:,si:si:sl) = A(:,si:si:sl)+ B(:,si:si:sl);
E(:,si) = 0.;

% --- define column indices for diagonals

Cm = [1:sl] - 1;
Dm = [1:sl];
Em = [1:sl] + 1;

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
    Lt = [ Lt, l+lo, l+lo, l+lo];
    Mt = [ Mt, Cm+mo, Dm+mo, Em+mo];
    Vt = [ Vt, C(b,:),D(b,:),E(b,:) ];
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

Hx = zeros(1,si);
Hy = zeros(1,si);

size(Hx)
size(Hy)

    for i=1:si
        l=i;
        Hx(i) = EVec(l,1);
        Hy(i) = EVec(l+sl,1);
    end


% --- calculate remaining field components

fprintf('calculating E field components ...\n');

Hz = zeros(1,si);
Ex = zeros(1,si);
Ey = zeros(1,si);

    for i=1:si-1
        Hz(i) = -1/keff*( (Hx(i+1)-Hx(i))/hx ); % = Hz / i
    end

Hz(si) = Hz(si-1);

emax = 0;

    for i=2:si
        l=i;
        Ex(i) = 1./(leps(l,1)*om) * (  keff*Hy(i) ); 
        Ey(i) = - 1./(leps(l,2)*om) * ( (Hz(i)-Hz(i-1))/hx + keff*Hx(i) ); 
        if abs(Ex(i)) > abs(emax) 
            emax = Ex(i);
        end
        if abs(Ey(i)) > abs(emax) 
            emax = Ey(i);
        end
    end


Ex(1) = Ex(2);
Ey(1) = Ey(2);

Ex = Ex/emax;
Ey = Ey/emax;
Hx = Hx/emax;
Hy = Hy/emax;

% --- plot fields

fprintf('plotting fields ...\n');

subplot(3,2,3);

hold on;
%contour(Ex',25);
title(sprintf('Ex component @ neff = %f',neff));
xlabel('x [a.u.]')
plot(Ex');
axis tight;
hold off;

subplot(3,2,4);

hold on;
%contour(Ey',25);
title(sprintf('Ey component @ neff = %f',neff));
xlabel('x [a.u.]')
plot(Ey');
hold off;

subplot(3,2,5);

hold on;
%contour(Hx',25);
title(sprintf('Hx component @ neff = %f',neff));
xlabel('x [a.u.]')
plot(Hx');
axis tight;
hold off;

subplot(3,2,6);
hold on;
title(sprintf('Hy component @ neff = %f',neff));
xlabel('x [a.u.]')
%contour(Hy',25);
plot(Hy');
axis tight;
hold off;

colormap jet(512);

fprintf('done.\n');






