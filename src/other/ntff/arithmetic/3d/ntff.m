function varargout = ntff(invlambda,freqnumber,theta,phi,eta0,varargin);
% varargout = ntff(invlambda,freqnumber,theta,phi,eta0,varargin)
% The function ntff computes the scattering cross sections of an object hit by a plane wave.
%
% The input is:
% - invlambda = inverse wavelength in computational units
% - freqnumber = the count for invlambda distinguishing files belonging to different invlambdas
% - theta = vector of spherical angle theta (inclination)
% - phi = vector of spherical angle phi (azimuth)
% - eta0 = free space impedance
% - varargin = variable input arguments (to keep backwards compatibility)
%   used to include the forward scattering direction and the field direction given as vector in spherical angles
%   [theta, phi, psi]. It also includes the distance between the NTFF plane and the 0 point. It can also
%   include the Courant factor DT
%
% When calling ntff the following files have to be present:
% dft%s_%i.set with %s element of {'+x','-x','+y','-y','+z','-z'}
%                   %i all numbers covered by freqnumber
% ../dft-ref_%i.set with %i same as above
%
% Sebastian Wuestner, 09.09.2009
% added extinction and absorption cross section on 06.11.2009

DT = 0.574; %Courant factor to collocate E and H in time (only minor difference if omitted)
if (size(varargin,2)==0)
   psi=-1;
elseif (size(varargin,2)==2)
   angles = varargin{1};
   theta_forw = angles(1);
   phi_forw = angles(2);
   psi = angles(3);
   tfsf_dist = varargin{2};
   n_theta_forw=find(theta==theta_forw);
   n_phi_forw=find(phi==phi_forw);
   if (size(n_theta_forw,1)==0 || size(n_phi_forw,1)==0) error('forward direction not probed'); end;
   l = sqrt(tfsf_dist(1)^2+tfsf_dist(3)^2)*cosd(atand(tfsf_dist(1)/tfsf_dist(3))-theta(n_theta_forw))-1;
   phase = exp(i*2*pi*invlambda*l);
elseif (size(varargin,2)==3)
   DT = varargin{1};
   angles = varargin{2};
   theta_forw = angles(1);
   phi_forw = angles(2);
   psi = angles(3);
   tfsf_dist = varargin{3};
   n_theta_forw=find(theta==theta_forw);
   n_phi_forw=find(phi==phi_forw);
   if (size(n_theta_forw,1)==0 || size(n_phi_forw,1)==0) error('forward direction not probed'); end;
   l = sqrt(tfsf_dist(1)^2+tfsf_dist(3)^2)*cosd(atand(tfsf_dist(1)/tfsf_dist(3))-theta(n_theta_forw))-1;
   phase = exp(i*2*pi*invlambda*l);
else
   error('Invalid number of variable input arguments; 0,2 or 3 input parameters.')
end;

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
for face = 1:6
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
         tmpexpR(1,:) = .5*tmpexpR(1,:); tmpexpR(size(tmpexpR,1),:) = .5*tmpexpR(size(tmpexpR,1),:);
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
Pinc = real(S(3))/2;

atheta = -i*kamp/(4*pi)*(Lphisum+eta0*Nthetasum);
aphi = i*kamp/(4*pi)*(Lthetasum-eta0*Nphisum);

%compute extinction cross section
if psi==0
   E = E * phase;
   Cext = imag(-4*pi/kamp*aphi(n_theta_forw,n_phi_forw)/E(2));% perpendicular (s)
elseif psi==90
   E = E * phase;
   Cext = imag(-4*pi/kamp*atheta(n_theta_forw,n_phi_forw)/E(2));% parallel (p)
elseif psi==-1
else
   error('Invalid psi-value! Extinction cross section cannot be evaluated.')
end;

%compute the RCS
RCS = 4*pi/(2*eta0*Pinc)*(atheta.*conj(atheta)+aphi.*conj(aphi));

switch nargout
   case 1
      %Differential scattering cross section
      varargout{1} = RCS;
   case 2
      %Differential scattering cross section
      varargout{1} = RCS;
      %Extinction cross section
      varargout{2} = Cext;
   otherwise
      error('Invalid number of output variables; 1 or 2 output variables.')
end;
