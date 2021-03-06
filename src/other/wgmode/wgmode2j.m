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

C = zeros(2,sl); % C -> i-1,j
D = zeros(2,sl); % D -> i,j
E = zeros(2,sl); % E -> i+1,j

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
        % --- block 2 (Hy-Hy)
        C(2,l) = leps(l,1)/leps(l,3)/(hx*hx);
        D(2,l) = leps(l,1)*om*om - leps(l,1)/(hx*hx)*(1./leps(l+1,3)+1./leps(l,3));
        E(2,l) = leps(l,1)/leps(l+1,3)/(hx*hx);    
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


keff = zeros(1,2);
neff = zeros(1,2);

lv = [1:sl];
Hc = zeros(2,si); 

for b = 1:2

Lt = [];
Mt = [];
Vt = [];

Lt = [ Lt, lv, lv, lv];
Mt = [ Mt, Cm, Dm, Em];
Vt = [ Vt, C(b,:),D(b,:),E(b,:) ];

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

S = sparse(L,M,V,sl,sl,length(V));


% --- run eigensolver

fprintf('running sparse solver ...\n');

opts.disp = 0;
[EVec,EVal] = eigs(S,[],1, betaguess,opts);

keff(1,b) = sqrt(EVal(1));

neff(1,b) = keff(1,b)/om;
fprintf(' -> keff = %f\n', keff(1,b));
fprintf(' -> neff = %f\n', neff(1,b));

    for i=1:si
        l=i;
        Hc(b,i) = EVec(l,1);
    end

end
    
Hx = Hc(1,:);
Hy = Hc(2,:);

% --- calculate remaining field components

fprintf('calculating E field components ...\n');

Hz = zeros(1,si);
Ex = zeros(1,si);
Ey = zeros(1,si);

    for i=1:si-1
        Hz(i) = -1/keff(1,1)*( (Hx(i+1)-Hx(i))/hx ); % = Hz / i
    end

Hz(si) = Hz(si-1);

emax = zeros(1,2);

    for i=2:si
        l=i;
        Ex(i) = 1./(leps(l,1)*om) * (  keff(1,2)*Hy(i) ); 
        Ey(i) = - 1./(leps(l,2)*om) * ( (Hz(i)-Hz(i-1))/hx + keff(1,1)*Hx(i) ); 
        if abs(Ex(i)) > abs(emax(1,2)) 
            emax(1,2) = Ex(i);
        end
        if abs(Ey(i)) > abs(emax(1,1)) 
            emax(1,1) = Ey(i);
        end
    end


Ex(1) = Ex(2);
Ey(1) = Ey(2);

Ex = Ex/emax(1,2);
Ey = Ey/emax(1,1);
Hx = Hx/emax(1,1);
Hy = Hy/emax(1,2);

% --- plot fields

fprintf('plotting fields ...\n');

subplot(3,2,3);

hold on;
%contour(Ex',25);
title(sprintf('Ex component @ neff = %f',neff(1,2)));
xlabel('x [a.u.]')
plot(Ex');
axis tight;
hold off;

subplot(3,2,4);

hold on;
%contour(Ey',25);
title(sprintf('Ey component @ neff = %f',neff(1,1)));
xlabel('x [a.u.]')
plot(Ey');
hold off;

subplot(3,2,5);

hold on;
%contour(Hx',25);
title(sprintf('Hx component @ neff = %f',neff(1,1)));
xlabel('x [a.u.]')
plot(Hx');
axis tight;
hold off;

subplot(3,2,6);
hold on;
title(sprintf('Hy component @ neff = %f',neff(1,2)));
xlabel('x [a.u.]')
%contour(Hy',25);
plot(Hy');
axis tight;
hold off;

colormap jet(512);

fprintf('done.\n');






