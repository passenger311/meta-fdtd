function RCS = ntff(invlambda,freqnumber,theta,phi,eta0);
% RCS = ntff(invlambda,freqnumber,theta,phi,eta0)
% The function ntff computes the differential scattering cross section of an object hit by a plane wave.
%  
% The input is:
% - invlambda = inverse wavelength in computational units
% - freqnumber = the count for invlambda distinguishing files belonging to different invlambdas
% - theta = vector of spherical angle theta (inclination)
% - phi = vector of spherical angle phi (azimuth)
% - eta0 = free space impedance
%
% When calling ntff the following files have to be present:
% dft%s_%i.set with %s element of {'+x','-x','+y','-y'}
%                   %i all numbers covered by freqnumber
% ../dft-ref_%i.set with %i same as above
%
% Sebastian Wuestner, 09.09.2009

DT = 0.7; %Courant factor to collocate E and H in time (only minor difference if omitted)

kamp = 2*pi*invlambda/eta0; % wavevector in medium with refractive index 1/eta0
max_theta = size(theta,2);
max_phi = size(phi,2);

%spherical projection vectors:
%array of unity vectors in r-direction er(theta,phi,component) where component = x,y,z
er = sind(theta)'*cosd(phi);
er(:,:,2) = sind(theta(:))*sind(phi(:))';
er(:,:,3) = cosd(theta(:))*ones(max_phi,1)';
%unity vector in theta-direction
etheta = cosd(theta)'*cosd(phi);
etheta(:,:,2) = cosd(theta(:))*sind(phi(:))';
etheta(:,:,3) = -sind(theta(:))*ones(max_phi,1)';
%unity vector in phi-direction
ephi = -ones(1,max_theta)'*sind(phi);
ephi(:,:,2) = ones(1,max_theta)'*cosd(phi);
ephi(:,:,3) = zeros(max_theta,max_phi);

%initialise arrays
Ltheta = zeros(max_theta,max_phi,6);
Lphi = zeros(max_theta,max_phi,6);
Ntheta = zeros(max_theta,max_phi,6);
Nphi = zeros(max_theta,max_phi,6);

%Sum over all faces of the NTFF-box
for face = 1:4
   %read in data from file
   [facerange, E, H] = readface(freqnumber,face);
   max_E = size(E,1);
   E = E * exp(i*pi*invlambda*DT);

   %surface currents via cross products
   J = zeros(max_E,3); M = zeros(max_E,3);
   if (face==1 | face ==2)
      n = (-1)^(face-1)*[ 1 0 0 ];
   elseif (face==3 | face ==4)
      n = (-1)^(face-1)*[ 0 1 0];
   elseif (face==5 | face ==6)
      n = (-1)^(face-1)*[ 0 0 1];
   end
   n = ones(max_E,1)*n;
   M = -cross(n,E,2);
   J = cross(n,H,2);

   %initialse arrays
   Jtheta = zeros(max_E,1); Mtheta = zeros(max_E,1);
   Jphi = zeros(max_E,1); Mphi = zeros(max_E,1);
   expR = zeros(max_E,1); tmpexpR = zeros(max_E,1);

   %sum over all theta-angles
   for n_theta = 1:max_theta
%      fprintf(1,'\ntheta = %f',theta(n_theta));
       %sum over all phi-angles
       for n_phi = 1:max_phi
%         fprintf(1,', phi = %f',phi(n_phi));
%         if mod(n_phi,5)==0 fprintf(1,'\n'); end

         ni = (facerange(1,1):facerange(1,3):facerange(1,2))*er(n_theta,n_phi,1);
         nj = (facerange(2,1):facerange(2,3):facerange(2,2))*er(n_theta,n_phi,2);
         nk = (facerange(3,1):facerange(3,3):facerange(3,2))*er(n_theta,n_phi,3);

         %exponential factors in integral
         if facerange(1,1)==facerange(1,2) %size(ni,2)==1
            tmpexpR = exp(i.*kamp.*(ni + nj'*ones(1,size(nk,2))+ones(size(nj,2),1)*nk));
         elseif facerange(2,1)==facerange(2,2) %size(nj,2)==1
            tmpexpR = exp(i.*kamp.*(ni'*ones(1,size(nk,2))+nj+ones(size(ni,2),1)*nk));
         elseif facerange(3,1)==facerange(3,2) %size(nk,2)==1
            tmpexpR = exp(i.*kamp.*(ni'*ones(1,size(nj,2))+ones(size(ni,2),1)*nj+nk));
         end
%         tmpexpR(1,:) = .5*tmpexpR(1,:); tmpexpR(size(tmpexpR,1),:) = .5*tmpexpR(size(tmpexpR,1),:);
         tmpexpR(:,1) = .5*tmpexpR(:,1); tmpexpR(:,size(tmpexpR,2)) = .5*tmpexpR(:,size(tmpexpR,2));
         expR=reshape(tmpexpR,max_E,1)*facerange(1,3)*facerange(2,3)*facerange(3,3);

         %projections of J and M on theta and phi unit vectors
         tmpetheta = squeeze(etheta(n_theta,n_phi,:));
         Jtheta = J*tmpetheta;
         Mtheta = M*tmpetheta;
         tmpephi = squeeze(ephi(n_theta,n_phi,:));
         Jphi = J*tmpephi;
         Mphi = M*tmpephi;

         %Ntheta etc by contracting Jtheta etc with expR
         Ntheta(n_theta,n_phi,face) = expR'*Jtheta;
         Ltheta(n_theta,n_phi,face) = expR'*Mtheta;
         Nphi(n_theta,n_phi,face) = expR'*Jphi;
         Lphi(n_theta,n_phi,face) = expR'*Mphi;
      end
   end
end
fprintf(1,'\n');

%take the sum over all faces of the cube
Nthetasum = squeeze(sum(Ntheta,3)); Nphisum = squeeze(sum(Nphi,3)); Lthetasum = squeeze(sum(Ltheta,3)); Lphisum = squeeze(sum(Lphi,3));

%normalisation
%reading reference file
clear E H;
[E,H] = readref(freqnumber);
E = E * exp(i*pi*invlambda*DT);

S = cross(E,conj(H));
Pinc = -real(S(3))/2;

%compute the RCS
RCS = 2*kamp/(eta0*Pinc).*(abs(Lphisum+eta0.*Nthetasum).*abs(Lphisum+eta0.*Nthetasum) + abs(Lthetasum-eta0.*Nphisum).*abs(Lthetasum-eta0.*Nphisum));
