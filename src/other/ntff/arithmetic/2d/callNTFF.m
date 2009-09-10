function callNTFF;

clear all

addpath /home/sw00100/Documents/MATLAB/NTFF_2d

fid = fopen('data.save','r');
if (fid==-1) error('file data.save does not exits'); end;
tline = fgetl(fid);
dx = sscanf(tline,'%e');
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

for freqnumber = 1:kmax
   invlambda = inv_lambda(freqnumber)
   RCS = ntff(invlambda,freqnumber,theta,phi,eta0);
   saveRCS('',freqnumber,theta,phi,dx*RCS);
   saveRCS2('',theta,phi,dx/invlambda,dx*RCS);
   RCS(:,1)=RCS(:,1)/2; RCS(:,size(phi,2))=RCS(:,size(phi,2))/2;
   fprintf(fid_Csca,'%e %e\n', dx/invlambda,dx/(2*pi)*(2*pi/(size(phi,2)-1))*sum(sind(theta)*RCS,2));
end
fclose(fid_Csca);

exit;
