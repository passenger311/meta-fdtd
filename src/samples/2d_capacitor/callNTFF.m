function callNTFF;

clear all

addpath /home/sw00100/Documents/MATLAB/NTFF_2d

fid = fopen('data.save','r');
if (fid==-1) error('file data.save does not exits'); end;
tline = fgetl(fid);
dx = sscanf(tline,'%e');
tline = fgetl(fid);
boxsize = sscanf(tline,'%e %e %e');
tline = fgetl(fid);
boxsize = sscanf(tline,'%e %e %e');
tline = fgetl(fid);
angles = sscanf(tline,'%e %e %e');
fclose(fid);

fid1 = fopen('invlambda.in', 'r');
if (fid1==-1) error('file invlambda.in does not exits'); end;
k=1;
while feof(fid1) ==0
   tline = fgetl(fid1);
   inv_lambda(k)=sscanf(tline,'%e');
   k=k+1;
end
fclose(fid1);
kmax=k-1;

theta = linspace(90,90,1);%0:90:180;
phi = linspace(-90,270,361);%0:180:360;
eta0 = 1;

fid_Csca = fopen('Csca_norm.dat','w');
fid_g = fopen('Csca_g.dat','w');
fid_g_geo = fopen('Csca_g_geo.dat','w');
fid_Cext = fopen('Cext.dat','w');
fid_Cabs = fopen('Cabs_norm.dat','w');

for freqnumber = 1:kmax
   invlambda = inv_lambda(freqnumber);
   [RCS,Cext] = ntff(invlambda,freqnumber,theta,phi,eta0,angles,boxsize./2+2);
   saveRCS('',freqnumber,theta,phi,RCS);
   RCS(:,1)=RCS(:,1)/2; RCS(:,size(phi,2))=RCS(:,size(phi,2))/2;
   Csca = dx/(2*pi)*(2*pi/(size(phi,2)-1))*sum(sind(theta)*RCS,2);
   fprintf(fid_Csca,'%e %e\n', dx/invlambda,Csca);
   fprintf(fid_g,'%e %e\n',dx/invlambda,sind(theta)*(RCS*cosd(phi-angles(2))')/sum(sind(theta)*RCS,2));
   fprintf(fid_g_geo,'%e %e\n',dx/invlambda,sind(theta)*(RCS*cosd(phi)')/sum(sind(theta)*RCS,2));
   fprintf(fid_Cext,'%e %e\n', dx/invlambda,dx*Cext);
   fprintf(fid_Cabs,'%e %e\n', dx/invlambda,dx*Cext-Csca);
end

fclose(fid_Csca);
fclose(fid_g);
fclose(fid_g_geo);
fclose(fid_Cext);
fclose(fid_Cabs);


%figure(1);
%surf(phi-(phi(2)-phi(1))/2,theta-(theta(2)-theta(1))/2,RCSsum)
%figure(1);
%surf(phi-(phi(2)-phi(1))/2,theta-(theta(2)-theta(1))/2,RCSsum)
%figure(1);
%surf(phi,theta,RCS);
%shading interp;
%figure(2);
%semilogy(theta,RCS);
exit;
