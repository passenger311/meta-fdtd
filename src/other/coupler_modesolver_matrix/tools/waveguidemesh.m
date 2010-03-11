function [x,y,xc,yc,nx,ny,eps] = waveguidemeshfull(eps_in,h,w,dx,dy);

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

jh = round(h/dy);
iw = round(w/dx);
jlayers = length(h);
ilayers = length(w);

nx = sum(iw);
ny = sum(jh)+1;

x = dx*(0:1:(nx))';
xc = (x(1:nx-1) + x(2:nx))/2;

y = (0:(ny-1))*dy;
yc = (1:(ny-1))*dy - dy/2;

eps = zeros(nx-1,ny-1);
ix=0;
for ii = 1:ilayers
  jy=0;
  for jj = 1:jlayers
    eps(ix+1:ix+iw(ii),jy+1:jy+jh(jj)) = eps_in(ii,jj);
    jy=jy+jh(jj);
  end
  ix=ix+iw(ii);
end

nx = length(xc);
ny = length(yc);

