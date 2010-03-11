function [x,y,xc,yc,nx,ny,eps] = waveguidemeshfull(n,h,rh,rw,gap,side,epsm,dx,dy);

% This function creates an index mesh for the finite-difference
% mode solver.  The function will accommodate a generalized three
% layer rib waveguide structure.  (Note: channel waveguides can
% also be treated by selecting the parameters appropriately.) 
% 
% USAGE:
% 
% [x,y,xc,yc,nx,ny,eps] = waveguidemeshfull(n,h,rh,rw,side,dx,dy)
%
% INPUT
%
% n - indices of refraction for layers in waveguide
% h - height of each layer in waveguide
% rh - height of waveguide feature
% rw - half-width of waveguide
% side - excess space to the right of waveguide
% dx - horizontal grid spacing
% dy - vertical grid spacing
% 
% OUTPUT
% 
% x,y - vectors specifying mesh coordinates
% xc,yc - vectors specifying grid-center coordinates
% nx,ny - size of index mesh
% eps - index mesh (n^2)
%
% AUTHOR:  Thomas E. Murphy (tem@umd.edu)

ih = round(h/dy);
irh = round (rh/dy);
irw = round (2*rw/dx);
igap = round (gap/dx);
iside = round (side/dx);
nlayers = length(h);

nx = irw+2*igap+2*iside+1;
ny = sum(ih)+1;

x = dx*(-(irw/2+igap+iside):1:(irw/2+igap+iside))';
xc = (x(1:nx-1) + x(2:nx))/2;

y = (0:(ny-1))*dy;
yc = (1:(ny-1))*dy - dy/2;

eps = zeros(nx-1,ny-1);

iy = 1;

for jj = 1:nlayers,
  for k = 1:ih(jj),
    eps(:,iy) = n(jj)^2*ones(nx-1,1);
    iy = iy+1;
  end
end

iy = ih(1)+ih(2)+1;%sum(ih)-ih(nlayers);
for k = 1:irh,
  eps(1:iside,iy) = (epsm(1)+i*epsm(2))*ones(iside,1);
  eps(irw+2*igap+iside+1:irw+2*igap+2*iside,iy) = (epsm(1)+i*epsm(2))*ones(iside,1);
  eps(iside+1:iside+igap,iy) = n(nlayers-1)^2*ones(igap,1);
  eps(irw+igap+iside+1:irw+2*igap+iside,iy) = n(nlayers-1)^2*ones(igap,1);
  iy = iy+1;
end

if (h(nlayers-1)~=0)
  ih(1)+ih(2)+ih(3)+1;
  for k = 1:ih(4),
    eps(1:iside,iy) = (epsm(1)+i*epsm(2))*ones(iside,1);
    eps(irw+2*igap+iside+1:irw+2*igap+2*iside,iy) = (epsm(1)+i*epsm(2))*ones(iside,1);
    iy = iy+1;
  end
  iy = sum(ih);
  for k = 1:ih(nlayers),
    eps(:,iy) = (epsm(1)+i*epsm(2))*ones(size(x,1)-1,1);
    eps(:,iy) = (epsm(1)+i*epsm(2))*ones(size(x,1)-1,1);
    iy = iy-1;
  end
end

nx = length(xc);
ny = length(yc);
