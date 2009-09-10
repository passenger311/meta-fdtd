function callNTFFbs;

clear all

addpath /home/sw00100/Documents/MATLAB/NTFF_final

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

theta = linspace(0,180,181);%0:90:180;
phi = linspace(0,90,2);%0:180:360;
eta0 = 1;

for freqnumber = 1:kmax
   invlambda = inv_lambda(freqnumber)
   RCS = ntffbs(invlambda,freqnumber,theta,phi,eta0);
   saveRCS('bs',freqnumber,theta,phi,RCS);
end

%figure(1);
%surf(phi-(phi(2)-phi(1))/2,theta-(theta(2)-theta(1))/2,RCSsum)
%figure(1);
%surf(phi-(phi(2)-phi(1))/2,theta-(theta(2)-theta(1))/2,RCSsum)
%figure(1);
%surf(phi,theta,RCS);
%shading interp;
%figure(2);
%semilogy(theta,RCS);

exit
