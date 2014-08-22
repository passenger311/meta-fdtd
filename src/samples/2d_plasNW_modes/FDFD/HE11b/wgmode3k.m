function [ Ex,Ey,Ez,Hx,Hy,Hz,range,neff] = wgmode3k( epsfile, lambdainv, hx, hy, neffguess, neff2, fig_on, cyl_rad, conv )
%  Detailed explanation goes here

% PML setup
pml=[1 1 0 1];	% 0:off, 1: on, xmin,xmax,ymin,ymax
d_pmlx=600;	% Thickness of PMLs in x direction
d_pmly=50;	% Thickness of PMLs in y direction
strength = 0.2;	% Strength of PMLs

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
tmp = textscan(fid, '%n %n %n');
fclose(fid);
leps = zeros(size(tmp{1},1),size(tmp,2));
for j = 1 : size(tmp,2)
   leps(:,j) = tmp{j};
end
clear tmp;

if ( si<1 || sj<1 )
 error('error')
end


fid = fopen(sprintf('%s_im_metal.in',epsfile),'r');

while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' ) 
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end

tmp = textscan(fid, '%n %n %n');
fclose(fid);
for j = 1 : size(tmp,2)
   leps(:,j) = leps(:,j)+i*tmp{j};
end
clear tmp;

leps_metal(:,1)=leps(:,1);
leps_metal(:,2)=leps(:,2);
leps_metal(:,3)=leps(:,3);

leps_metal(sl+1:sl+si,1) = leps_metal(sl,1); % add a row
leps_metal(sl+1:sl+si,2) = leps_metal(sl,2); % add a row
leps_metal(sl+1:sl+si,3) = leps_metal(sl,3); % add a row 

fid = fopen(sprintf('%s_im_gain.in',epsfile),'r');

while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' ) 
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end

tmp = textscan(fid, '%n %n %n');
fclose(fid);
for j = 1 : size(tmp,2)
   leps_gain(:,j) = i*tmp{j};
   leps(:,j) = leps(:,j)+leps_gain(:,j);
end
clear tmp;

leps_gain(sl+1:sl+si,1) = leps_gain(sl,1); % add a row
leps_gain(sl+1:sl+si,2) = leps_gain(sl,2); % add a row
leps_gain(sl+1:sl+si,3) = leps_gain(sl,3); % add a row 

% Set up PMLs

scalx0=zeros(si,sj); scalxh0=scalx0; scaly0=scalx0; scalyh0=scalx0;
if (pml(1)==1);
  for k=1:d_pmlx;
    scalx0(k,:)=i*strength*(d_pmlx-(k-1)).^2*ones(1,sj)/d_pmlx^2;
    scalxh0(k,:)=i*strength*(d_pmlx-(k-0.5)).^2*ones(1,sj)/d_pmlx^2;
  end
end
if (pml(2)==1);
  for k=1:d_pmlx;
    scalx0(end-(k-1),:)=i*strength*(d_pmlx-(k-1)).^2*ones(1,sj)/d_pmlx^2;
    scalxh0(end-(k-1),:)=i*strength*(d_pmlx-(k-1.5)).^2*ones(1,sj)/d_pmlx^2;
  end
end
if (pml(3)==1);
  for k=1:d_pmly;
    scaly0(:,k)=i*strength*(d_pmly-(k-1)).^2*ones(si,1)/d_pmly^2;
    scalyh0(:,k)=i*strength*(d_pmly-(k-0.5)).^2*ones(si,1)/d_pmly^2;
  end
end
if (pml(4)==1);
  for k=1:d_pmly;
    scaly0(:,end-(k-1))=i*strength*(d_pmly-(k-1)).^2*ones(si,1)/d_pmly^2;
    scalyh0(:,end-(k-1))=i*strength*(d_pmly-(k-1.5)).^2*ones(si,1)/d_pmly^2;
  end
end

scalxh=reshape(scalxh0,si*sj,1);
scalx=reshape(scalx0,si*sj,1);
scalyh=reshape(scalyh0,si*sj,1);
scaly=reshape(scaly0,si*sj,1);



leps(:,1)=leps(:,1)+scalxh+scaly;
leps(:,2)=leps(:,2)+scalx+scalyh;
leps(:,3)=leps(:,3)+scalxh+scalyh;

leps(sl+1:sl+si,1) = leps(sl,1); % add a row
leps(sl+1:sl+si,2) = leps(sl,2); % add a row
leps(sl+1:sl+si,3) = leps(sl,3); % add a row 


fprintf('plotting epsilon ...\n');

for j=1:3;
  eps2(:,:,j)=conj(reshape(leps(1:sl,j),si,sj)');
  eps_metal(:,:,j)=conj(reshape(leps_metal(1:sl,j),si,sj)');
  eps_gain(:,:,j)=conj(reshape(leps_gain(1:sl,j),si,sj)');
end
nmax = max(max(abs(sqrt(leps))));

fprintf(' -> nmax = %f\n', nmax);


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

% --- block 1 (Hx-Hx)
A(1,:) = leps(1:sl,2)./leps(1:sl,3)./(hy.*hy);
C(1,:) = 1./(hx.*hx);
D(1,:) = -1./(hx.*hx) - 1./(hx.*hx) + leps(1:sl,2).*om.*om - leps(1:sl,2)./(hy.*hy).*(1./leps(si+1:sl+si,3)+1./leps(1:sl,3));
E(1,:) = 1./(hx.*hx);
G(1,:) = leps(1:sl,2)./leps(si+1:sl+si,3)./(hy.*hy); 
% --- block 2 (Hx-Hy)
C(2,:) = 1./(hx.*hy).*(1-leps(1:sl,2)./leps(1:sl,3));
D(2,:) = -C(2,:);
F(2,:) = -1./(hx.*hy).*(1-leps(1:sl,2)./leps(si+1:sl+si,3));
G(2,:) = -F(2,:);
% --- block 3 (Hy-Hy)
A(3,:) = 1./(hy.*hy);
C(3,:) = leps(1:sl,1)./leps(1:sl,3)./(hx.*hx);
D(3,:) = -1./(hy.*hy) -1./(hy.*hy) + leps(1:sl,1)*om*om - leps(1:sl,1)./(hx.*hx).*(1./leps(2:sl+1,3)+1./leps(1:sl,3));
E(3,:) = leps(1:sl,1)./leps(2:sl+1,3)./(hx.*hx);
G(3,:) = 1./(hy.*hy);
% --- block 4 (Hy-Hx)
A(4,:) = 1./(hx.*hy).*(1-leps(1:sl,1)./leps(1:sl,3));
B(4,:) = -1./(hx.*hy).*(1-leps(1:sl,1)./leps(2:sl+1,3));
D(4,:) = -A(4,:);
E(4,:) = -B(4,:);       


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


fprintf('wiping out zero values ...\n');

L=Lt(Vt~=0);
M=Mt(Vt~=0);
V=Vt(Vt~=0);

clear A B C D E F G Lt Mt Vt

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

if(exist('mode_neff.dat')~=2)
  f_out=fopen('mode_neff.dat','w');
  fprintf(f_out,'#si\t#lambda\tRe neff\tIm neff\tProp. length (micron)\n');
  fclose(f_out);
end
f_out=fopen('mode_neff.dat','a');
fprintf(f_out,'%4.0f\t%3.0f\t%3.2f\t%5.4f\t%5.4f\t%f\n',si,2*cyl_rad,conv/lambdainv,real(neff),imag(neff),conv/lambdainv/2/pi/imag(neff)/1000);
fclose(f_out);

f_out=fopen('neff.lua','w');
fprintf(f_out,'nrefr = %f\n',neff);
fprintf(f_out,'nrefr_im = %f\n',imag(neff));
fclose(f_out);

Hx = conj(reshape(EVec(1:sl,1),si,sj)');
Hy = conj(reshape(EVec(sl+1:2*sl,1),si,sj)');


% --- calculate remaining field components

fprintf('calculating E field components ...\n');

Hz(1:sj-1,1:si-1) = -i./keff.*( (Hx(1:sj-1,2:si)-Hx(1:sj-1,1:si-1))./hx + (Hy(2:sj,1:si-1)-Hy(1:sj-1,1:si-1))./hy); % = Hz / ii

Hz(sj,:) = Hz(sj-1,:);
Hz(:,si) = Hz(:,si-1);

Ex(2:sj,2:si) = 1./(eps2(2:sj,2:si,1).*om) .* ( -i.*(Hz(2:sj,2:si)-Hz(1:sj-1,2:si))./hy + keff.*Hy(2:sj,2:si) );
Ey(2:sj,2:si) = - 1./(eps2(2:sj,2:si,2).*om) .* ( -i.*(Hz(2:sj,2:si)-Hz(2:sj,1:si-1))./hx + keff.*Hx(2:sj,2:si) );
Ez(2:sj,2:si) = i./(eps2(2:sj,2:si,3).*om) .* ( (Hx(2:sj,2:si)-Hx(1:sj-1,2:si))./hy - (Hy(2:sj,2:si)-Hy(2:sj,1:si-1))./hx );
emax = max(max(max(Ex)),max(max(Ey)));

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

fid = fopen('geo_nanowire.in','r');
while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' )
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end
tmp = textscan(fid, '%n %n %n');
fclose(fid);
leps = zeros(size(tmp{1},1),size(tmp,2));
for j = 1 : size(tmp,2)
   leps_nw(:,j) = tmp{j};
end
clear tmp;
leps_nw(sl+1:sl+si,1) = leps_nw(sl,1); % add a row
leps_nw(sl+1:sl+si,2) = leps_nw(sl,2); % add a row
leps_nw(sl+1:sl+si,3) = leps_nw(sl,3); % add a row 

for j=1:3;
  eps_nw(:,:,j)=conj(reshape(leps_nw(1:sl,j),si,sj)');
end

fid = fopen('geo_energy.in','r');
while feof(fid) == 0
    tline = fgetl(fid);
    if ( tline(1) =='!' || tline(2)=='!' || tline(1) == '(' )
        continue;
    end
    range = sscanf(tline,'%i %i %i %i %i %i %i %i %i')';
    break;
end
tmp = textscan(fid, '%n %n %n');
fclose(fid);
leps = zeros(size(tmp{1},1),size(tmp,2));
for j = 1 : size(tmp,2)
   leps_en(:,j) = tmp{j};
end
clear tmp;
leps_en(sl+1:sl+si,1) = leps_en(sl,1); % add a row
leps_en(sl+1:sl+si,2) = leps_en(sl,2); % add a row
leps_en(sl+1:sl+si,3) = leps_en(sl,3); % add a row 

for j=1:3;
  eps_en(:,:,j)=conj(reshape(leps_en(1:sl,j),si,sj)');
end


% 2 times electric energy in gain area
En_g=2*sum(sum(eps_nw(:,:,1).*abs(Ex(:,:)).^2+eps_nw(:,:,2).*abs(Ey(:,:)).^2+eps_nw(:,:,3).*abs(Ez(:,:)).^2,1),2);
% Total energy
En_tot=sum(sum(eps_en(:,:,1).*abs(Ex(:,:)).^2+eps_en(:,:,2).*abs(Ey(:,:)).^2+eps_en(:,:,3).*abs(Ez(:,:)).^2,1),2)+sum(sum(abs(Hx(:,:)).^2+abs(Hy(:,:)).^2+abs(Hz(:,:)).^2,1),2);
% Confinement factor
Gamma = En_g/En_tot;
A_g = sum(sum((eps_nw(:,:,1)/max(max(max(eps_nw))))>0,2),1)*conv^2;
% Effective area defined using confinement factor
A_eff0 = A_g/Gamma;
fprintf('\nConfinement factor and effective mode area: %7.4f %3.1e\n',Gamma,A_eff0)

% Energy squared
En2_tot=sum(sum((eps_en(:,:,1).*abs(Ex(:,:)).^2+eps_en(:,:,2).*abs(Ey(:,:)).^2+eps_en(:,:,3).*abs(Ez(:,:)).^2).^2,1),2)+sum(sum(abs(Hx(:,:)).^2+abs(Hy(:,:)).^2+abs(Hz(:,:)).^2,1),2);
% Maximum local energy
max_En = max(max(eps_en(:,:,1).*abs(Ex(:,:)).^2+eps_en(:,:,2).*abs(Ey(:,:)).^2+eps_en(:,:,3).*abs(Ez(:,:)).^2));
% Effective area defined using energy squared
A_eff_En2 = En_tot^2/En2_tot*conv^2;
% Effective area defined using maximum local energy
A_eff_Emax = En_tot/max_En*conv^2;

fprintf('Alternative measures of the effective mode area: En^2: %3.1e, max_En: %3.1e\n',A_eff_En2,A_eff_Emax)

% Rate defined using average energy velocity

inside_x=1:size(eps2,2);
inside_x=d_pmlx+1:size(eps2,2)-d_pmlx-1;
inside_y=1:size(eps2,1);
inside_y=1:size(eps2,1)-d_pmly-1;
En_tot_inside=eps_en(:,:,1).*abs(Ex(:,:)).^2+eps_en(:,:,2).*abs(Ey(:,:)).^2+eps_en(:,:,3).*abs(Ez(:,:)).^2+abs(Hx(:,:)).^2+abs(Hy(:,:)).^2+abs(Hz(:,:)).^2;
En_tot_inside=sum(sum(En_tot_inside(inside_y,inside_x)));
gamma_inside=om*(2.9979e8/conv/1e-9)/1e0*2*0.5*imag((eps2(:,:,1).*abs(Ex(:,:)).^2+eps2(:,:,2).*abs(Ey(:,:)).^2+eps2(:,:,3).*abs(Ez(:,:)).^2)/(0.5*En_tot_inside));
gamma_inside=sum(sum(gamma_inside(inside_y,inside_x)));
Sz=real(Ex(:,:).*conj(Hy(:,:))-Ey(:,:).*conj(Hx(:,:)));
Sz=sum(sum(Sz(inside_y,inside_x)));
n_Eav_inside=(0.5*En_tot_inside)./Sz;
Sx=real(Ey(:,:).*conj(Hz(:,:))-Ez(:,:).*conj(Hy(:,:)));
Sx=-sum(Sx(inside_y,inside_x(1)))+sum(Sx(inside_y,inside_x(end)));
Sy=real(Ez(:,:).*conj(Hx(:,:))-Ex(:,:).*conj(Hz(:,:)));
Sy=-sum(Sy(inside_y(1),inside_x))+sum(Sy(inside_y(end),inside_x));
alpha_Sx=-Sx/Sz/conv*1e7;
alpha_Sy=-Sy/Sz/conv*1e7;
alpha_inside= gamma_inside*n_Eav_inside/2.9979e8/100;

gamma=om*(2.9979e8/conv/1e-9)/1e0*2*0.5*imag(sum(sum(eps2(:,:,1).*abs(Ex(:,:)).^2+eps2(:,:,2).*abs(Ey(:,:)).^2+eps2(:,:,3).*abs(Ez(:,:)).^2,1),2)/(0.5*En_tot));
gamma_metal=om*(2.9979e8/conv/1e-9)/1e0*2*0.5*imag(sum(sum(eps_metal(:,:,1).*abs(Ex(:,:)).^2+eps_metal(:,:,2).*abs(Ey(:,:)).^2+eps_metal(:,:,3).*abs(Ez(:,:)).^2,1),2)/(0.5*En_tot));
gamma_gain=om*(2.9979e8/conv/1e-9)/1e0*2*0.5*imag(sum(sum(eps_gain(:,:,1).*abs(Ex(:,:)).^2+eps_gain(:,:,2).*abs(Ey(:,:)).^2+eps_gain(:,:,3).*abs(Ez(:,:)).^2,1),2)/(0.5*En_tot));
% Average energy velocity I/W (different from group velocity) 
n_Eav=(0.5*En_tot)/real(sum(sum(Ex(:,:).*conj(Hy(:,:))-Ey(:,:).*conj(Hx(:,:)),1),2));
% Loss coefficient from fields and flux
alpha=gamma*n_Eav/2.9979e8/100;
alpha_metal=gamma_metal*n_Eav/2.9979e8/100;
alpha_gain=gamma_gain*n_Eav/2.9979e8/100;
% Loss coefficient from effective index
alpha2=2*om*(1/conv/1e-9)*imag(neff)/100;
gamma = alpha2/n_Eav*2.9979e8*100;

if(exist('mode_full.dat')~=2)
  f_out=fopen('mode_full.dat','w');
  fprintf(f_out,'#si\t#Radius\tlambda\tRe neff\tIm neff\talpha\talph in\talph Sx\talph Sy\talph me\talph ga\tConf\tA_0\tA_En2\tA_Emax\tn_Eav\tgamma\talpha_g_bulk\n');
  fclose(f_out);
%  f_out=fopen('mode_alt_area.dat','w');
%  fprintf(f_out,'#lambda\tRe neff \tA_0\tA_En2\tA_Emax\n');
%  fclose(f_out);
%  f_out=fopen('mode_gamma_area.dat','w');
%  fprintf(f_out,'#lambda\tRe(neff)\tConf\tA_0\n');
%  fclose(f_out);
end
%  f_out=fopen('mode_gamma_area.dat','a');
%  fprintf(f_out,'%3.2f\t%4.3f\t%f\t%6.0f\n',conv/lambdainv,real(neff),Gamma,A_eff0);
%  fclose(f_out);
%  f_out=fopen('mode_alt_area.dat','a');
%  fprintf(f_out,'%3.2f\t%4.3f\t%8.0f\t%6.0f\t%6.0f\n',conv/lambdainv,real(neff),A_eff0,A_eff_En2,A_eff_Emax);
%  fclose(f_out);
f_out=fopen('mode_full.dat','a');
fprintf(f_out,'%4.0f\t%3.0f\t%3.2f\t%5.4f\t%5.4f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%4.3f\t%6.0f\t%6.0f\t%6.0f\t%6.2f\t%6.2f\t%6.0f\n',si,2*cyl_rad,conv/lambdainv,real(neff),imag(neff),alpha,alpha_inside,alpha_Sx,alpha_Sy,alpha_metal,alpha_gain,Gamma,A_eff0,A_eff_En2,A_eff_Emax,n_Eav,gamma/1e12,alpha2*2.3/n_Eav/Gamma);
fclose(f_out);


% --- plot fields

if (fig_on==1);
fprintf('plotting fields ...\n');
figure('visible','off')
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
pcolor(real(eps2(:,:,1)));
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

filename=sprintf('mode_%4.0f_%3.0f_%3.0f_%5.3f+i%5.3f.png',si,2*cyl_rad,conv/lambdainv,real(neff),imag(neff));
print(gcf, '-dpng', '-r150', filename);
close(figure);
end

fprintf('done.\n');






