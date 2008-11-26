function [phix,phiy,phiz,psix,psiy,psiz,neff,blob] = wgmodes (fileName,lambda, guess, nmodes, boundary);


% This function computes the two transverse magnetic field
% components of a dielectric waveguide, using the finite
% difference method.  For details about the method, please
% consult:  
%
% A. B. Fallahkhair, K. S. Li, and T. E. Murphy, "Vector
% Finite Difference Modesolver for Anisotropic Dielectric
% Waveguides", submitted to J. Lightwave. Technol., 2007.
%
% USAGE:
% 
% [hx,hy,neff] = wgmodes(lambda, guess, nmodes, dx, dy, ...
%                        eps,boundary);
% [hx,hy,neff] = wgmodes(lambda, guess, nmodes, dx, dy, ...
%                        epsxx, epsyy, epszz, boundary);
% [hx,hy,neff] = wgmodes(lambda, guess, nmodes, dx, dy, ...
%                        epsxx, epsxy, epsyx, epsyy, epszz, boundary);
% 
% INPUT:
% 
% lambda - optical wavelength
% guess - scalar shift to apply when calculating the eigenvalues.
%     This routine will return the eigenpairs which have an
%     effective index closest to this guess
% nmodes - the number of modes to calculate
% dx - horizontal grid spacing (vector or scalar)
% dy - vertical grid spacing (vector or scalar)
% eps - index mesh (isotropic materials)  OR:
% epsxx, epsxy, epsyx, epsyy, epszz - index mesh (anisotropic)
% boundary - 4 letter string specifying boundary conditions to be
% applied at the edges of the computation window.  
%   boundary(1) = North boundary condition
%   boundary(2) = South boundary condition
%   boundary(3) = East boundary condition
%   boundary(4) = West boundary condition
% The following boundary conditions are supported: 
%   'A' - Hx is antisymmetric, Hy is symmetric.
%   'S' - Hx is symmetric and, Hy is antisymmetric.
%   '0' - Hx and Hy are zero immediately outside of the
%         boundary. 
% 
% OUTPUT:
% 
% hx - three-dimensional vector containing Hx for each
%      calculated mode 
% hy - three-dimensional vector containing Hy for each
%      calculated mode (e.g.: hy(:,k) = two dimensional Hy
%      matrix for the k-th mode 
% neff - vector of modal effective indices
%
% NOTES:
%
% 1) The units are arbitrary, but they must be self-consistent
% (e.g., if lambda is in um, then dx and dy should also be in
% um.
%
% 2) Unlike the E-field modesolvers, this method calculates
% the transverse MAGNETIC field components Hx and Hy.  Also,
% it calculates the components at the edges (vertices) of
% each cell, rather than in the center of each cell.  As a
% result, if size(eps) = [n,m], then the output eigenvectors
% will be have a size of [n+1,m+1].
%
% 3) This version of the modesolver can optionally support
% non-uniform grid sizes.  To use this feature, you may let dx
% and/or dy be vectors instead of scalars.
%
% 4) The modesolver can consider anisotropic materials, provided
% the permittivity of all constituent materials can be
% expressed in one of the following forms:   
%
%  [eps  0   0 ]  [epsxx   0     0  ]  [epsxx epsxy   0  ]
%  [ 0  eps  0 ]  [  0   epsyy   0  ]  [epsyx epsyy   0  ]
%  [ 0   0  eps]  [  0     0   epszz]  [  0     0   epszz]
%
% The program will decide which form is appropriate based upon
% the number of input arguments supplied.
%
% 5) Perfectly matched boundary layers can be accomodated by
% using the complex coordinate stretching technique at the
% edges of the computation window.  (stretchmesh.m can be used
% for complex or real-coordinate stretching.)
%
% AUTHORS:  Thomas E. Murphy (tem@umd.edu)
%           Arman B. Fallahkhair (a.b.fallah@gmail.com)
%           Kai Sum Li (ksl3@njit.edu)

%open file
fid = fopen(fileName, 'r');
if (fid==-1)
  fprintf('\nerror: cannot open %s.\n', fileName);
end
fprintf(1,'reading file %s...\n', fileName);

while feof(fid) == 0
   tline = fgetl(fid);
   if (length(tline)==0)
     continue;
   end
   if (tline(1)=='(')|(tline(2)=='!')
     continue;
   end
   [blob] = sscanf(tline,'%i %i %i %i %i %i %i %i %i');
   break;
end;

dx = blob(3);
nrow = ceil((blob(2)-blob(1))/dx+.5);
dy = blob(6);
ncolumn = ceil((blob(5)-blob(4))/dy+.5);
%read in file
epsxx=zeros(nrow,ncolumn);epsyy=zeros(nrow,ncolumn);
epszz=zeros(nrow,ncolumn);
row=1;column=1;

while feof(fid) == 0
   tline = fgetl(fid);
   if (length(tline)==0)
     continue;
   end
   if (tline(1)==')')
     continue;
   end
   tmp = sscanf(tline,'%e %e %e');
   epsxx(row,column)=tmp(1);
   epsyy(row,column)=tmp(2);
   epszz(row,column)=tmp(3);
   row = row+1;
   if(row>nrow);
     column=column+1;
     row=1;
   end;
end
fclose(fid);




fprintf (1,'solving for eigenmodes...'); t = cputime;



%if (nargin == 11)
%  epsxx = varargin{1};
%  epsxy = varargin{2};
%  epsyx = varargin{3};
%  epsyy = varargin{4};
%  epszz = varargin{5};
%  boundary = varargin{6};
%  if (nargin == 8)
%    epsxx = varargin{1};
%    epsxy = zeros(size(epsxx));
%    epsyx = zeros(size(epsxx));
%    epsyy = varargin{2};
%    epszz = varargin{3};
%    boundary = varargin{4};
%  elseif (nargin == 6)
%    epsxx = varargin{1};
%    epsxy = zeros(size(epsxx));
%    epsyx = zeros(size(epsxx));
%    epsyy = epsxx;
%  epszz = epsxx;
%  boundary = varargin{2};
%else
%  error('Incorrect number of input arguments.\n');
%end

[nx,ny] = size(epsxx);
%nx = nx + 1;
%ny = ny + 1;

% now we pad eps on all sides by one grid point
%epsxx = [epsxx(:,1),epsxx,epsxx(:,ny-1)];
%epsxx = [epsxx(1,1:ny+1);epsxx;epsxx(nx-1,1:ny+1)];

%epsyy = [epsyy(:,1),epsyy,epsyy(:,ny-1)];
%epsyy = [epsyy(1,1:ny+1);epsyy;epsyy(nx-1,1:ny+1)];

%epsxy = [epsxy(:,1),epsxy,epsxy(:,ny-1)];
%epsxy = [epsxy(1,1:ny+1);epsxy;epsxy(nx-1,1:ny+1)];

%epsyx = [epsyx(:,1),epsyx,epsyx(:,nysize(epsxx)-1)];
%epsyx = [epsyx(1,1:ny+1);epsyx;epsyx(nx-1,1:ny+1)];

%epszz = [epszz(:,1),epszz,epszz(:,ny-1)];
%epszz = [epszz(1,1:ny+1);epszz;epszz(nx-1,1:ny+1)];
epszz = [epszz,epszz(:,ny)];
epszz = [epszz;epszz(nx,1:ny+1)];

k = 2*pi/lambda;  % free-space wavevector

%if isscalar(dx)
%  dx = dx*ones(nx,1);             % uniform grid
%else
%  dx = dx(:);                       % convert to column vector
%  dx = [dx(1);dx;dx(length(dx))];   % pad dx on top and bottom
%end

%if isscalar(dy)
%  dy = dy*ones(1,ny);             % uniform grid
%else
%  dy = dy(:);                       % convert to column vector
%  dy = [dy(1);dy;dy(length(dy))]';  % pad dy on top and bottom
%end

% distance to neighboring points to north south east and west,
% relative to point under consideration (P), as shown below.

%n = ones(1,nx*ny);      n(:) = ones(nx,1)*dy(2:ny+1);
%s = ones(1,nx*ny);      s(:) = ones(nx,1)*dy(1:ny);
%e = ones(1,nx*ny);      e(:) = dx(2:nx+1)*ones(1,ny);
%w = ones(1,nx*ny);      w(:) = dx(1:nx)*ones(1,ny);

% epsilon tensor elements in regions 1,2,3,4, relative to the
% mesh point under consideration (P), as shown below.
%
%                 NW------N------NE
%                 |       |       |
%                 |   1   n   4   |
%                 |       |       |
%                 W---w---P---e---E
%                 |       |       |
%                 |   2   s   3   |
%                 |       |       |
%                 SW------S------SE

%exx1 = ones(1,nx*ny);   exx1(:) = epsxx(1:nx,2:ny+1);
%exx2 = ones(1,nx*ny);   exx2(:) = epsxx(1:nx,1:ny);
%exx3 = ones(1,nx*ny);   exx3(:) = epsxx(2:nx+1,1:ny);
%exx4 = ones(1,nx*ny);   exx4(:) = epsxx(2:nx+1,2:ny+1);
exx = ones(1,nx*ny);   exx(:) = epsxx(1:nx,1:ny);


%eyy1 = ones(1,nx*ny);   eyy1(:) = epsyy(1:nx,2:ny+1);
%eyy2 = ones(1,nx*ny);   eyy2(:) = epsyy(1:nx,1:ny);
%eyy3 = ones(1,nx*ny);   eyy3(:) = epsyy(2:nx+1,1:ny);
%eyy4 = ones(1,nx*ny);   eyy4(:) = epsyy(2:nx+1,2:ny+1);
eyy = ones(1,nx*ny);   eyy(:) = epsyy(1:nx,1:ny);

%exy1 = ones(1,nx*ny);   exy1(:) = epsxy(1:nx,2:ny+1);
%exy2 = ones(1,nx*ny);   exy2(:) = epsxy(1:nx,1:ny);
%exy3 = ones(1,nx*ny);   exy3(:) = epsxy(2:nx+1,1:ny);
%exy4 = ones(1,nx*ny);   exy4(:) = epsxy(2:nx+1,2:ny+1);

%eyx1 = ones(1,nx*ny);   eyx1(:) = epsyx(1:nx,2:ny+1);
%eyx2 = ones(1,nx*ny);   eyx2(:) = epsyx(1:nx,1:ny);
%eyx3 = ones(1,nx*ny);   eyx3(:) = epsyx(2:nx+1,1:ny);
%eyx4 = ones(1,nx*ny);   eyx4(:) = epsyx(2:nx+1,2:ny+1);

%ezz1 = ones(1,nx*ny);   ezz1(:) = epszz(1:nx,2:ny+1);
%ezz2 = ones(1,nx*ny);   ezz2(:) = epszz(1:nx,1:ny);
%ezz3 = ones(1,nx*ny);   ezz3(:) = epszz(2:nx+1,1:ny);
%ezz4 = ones(1,nx*ny);   ezz4(:) = epszz(2:nx+1,2:ny+1);
ezzy = ones(1,nx*ny);   ezzy(:) = epszz(1:nx,2:ny+1);
ezzx = ones(1,nx*ny);   ezzx(:) = epszz(2:nx+1,1:ny);
ezz = ones(1,nx*ny);   ezz(:) = epszz(1:nx,1:ny);

%ns21 = n.*eyy2+s.*eyy1;
%ns34 = n.*eyy3+s.*eyy4;
%ew14 = e.*exx1+w.*exx4;
%ew23 = e.*exx2+w.*exx3;


%axxn = ((2*eyy4.*e-eyx4.*n).*(eyy3./ezz4)./ns34 + ...
%        (2*eyy1.*w+eyx1.*n).*(eyy2./ezz1)./ns21)./(n.*(e+w));
axxn = (eyy./ezzy)./dy^2;


%axxs = ((2*eyy3.*e+eyx3.*s).*(eyy4./ezz3)./ns34 + ...
%        (2*eyy2.*w-eyx2.*s).*(eyy1./ezz2)./ns21)./(s.*(e+w));
axxs = eyy./ezz./dy^2;


%ayye = (2.*n.*exx4 - e.*exy4).*exx1./ezz4./e./ew14./(n+s) + ...
%       (2.*s.*exx3 + e.*exy3).*exx2./ezz3./e./ew23./(n+s);
ayyn = ones(1,nx*ny)./dy^2;

%ayyw = (2.*exx1.*n + exy1.*w).*exx4./ezz1./w./ew14./(n+s) + ...
%       (2.*exx2.*s - exy2.*w).*exx3./ezz2./w./ew23./(n+s);
ayys = ones(1,nx*ny)./dy^2;

%axxe = 2./(e.*(e+w)) + ...
%       (eyy4.*eyx3./ezz3 - eyy3.*eyx4./ezz4)./(e+w)./ns34;
axxe = ones(1,nx*ny)./dx^2;

%axxw = 2./(w.*(e+w)) + ...
%       (eyy2.*eyx1./ezz1 - eyy1.*eyx2./ezz2)./(e+w)./ns21;
axxw = ones(1,nx*ny)./dx^2;

%ayyn = 2./(n.*(n+s)) + ...
%       (exx4.*exy1./ezz1 - exx1.*exy4./ezz4)./(n+s)./ew14;
ayye = exx./ezzx./dx^2;

%ayys = 2./(s.*(n+s)) + ...
%       (exx2.*exy3./ezz3 - exx3.*exy2./ezz2)./(n+s)./ew23;
ayyw = exx./ezz./dx^2;

%axxne = +eyx4.*eyy3./ezz4./(e+w)./ns34;
%axxse = -eyx3.*eyy4./ezz3./(e+w)./ns34;
%axxnw = -eyx1.*eyy2./ezz1./(e+w)./ns21;
%axxsw = +eyx2.*eyy1./ezz2./(e+w)./ns21;
%axxne = zeros(1,nx*ny); axxse = zeros(1,nx*ny);
%axxnw = zeros(1,nx*ny); axxsw = zeros(1,nx*ny);

%ayyne = +exy4.*exx1./ezz4./(n+s)./ew14;
%ayyse = -exy3.*exx2./ezz3./(n+s)./ew23;
%ayynw = -exy1.*exx4./ezz1./(n+s)./ew14;
%ayysw = +exy2.*exx3./ezz2./(n+s)./ew23;
%ayyne = zeros(1,nx*ny); ayyse = zeros(1,nx*ny);
%ayynw = zeros(1,nx*ny); ayysw = zeros(1,nx*ny);

%axxp = - axxn - axxs - axxe - axxw - axxne - axxse - axxnw - axxsw ...
%       + k^2*(n+s).*(eyy4.*eyy3.*e./ns34 + eyy1.*eyy2.*w./ns21)./(e+w);
axxp = -2.0.*ones(1,nx*ny)./dx^2 + eyy.*k^2 - eyy./dy^2.*(1./ezzy+1./ezz);

%ayyp = - ayyn - ayys - ayye - ayyw - ayyne - ayyse - ayynw - ayysw ...
%       + k^2*(e+w).*(exx1.*exx4.*n./ew14 + exx2.*exx3.*s./ew23)./(n+s);
ayyp = -2.0.*ones(1,nx*ny)./dy^2 + exx.*k^2 - exx./dx^2.*(1./ezzx+1./ezz);

%axyn = (eyy3.*eyy4./ezz4./ns34 - ...
%        eyy2.*eyy1./ezz1./ns21 + ...
%        s.*(eyy2.*eyy4 - eyy1.*eyy3)./ns21./ns34)./(e+w);
axyn = (ones(1,nx*ny)-eyy./ezzy)./dx./dy;

%axys = (eyy1.*eyy2./ezz2./ns21 - ...
%        eyy4.*eyy3./ezz3./ns34 + ...
%        n.*(eyy2.*eyy4 - eyy1.*eyy3)./ns21./ns34)./(e+w);
axys = zeros(1,nx*ny);

%ayxe = (exx1.*exx4./ezz4./ew14 - ...
%        exx2.*exx3./ezz3./ew23 + ...
%        w.*(exx2.*exx4 - exx1.*exx3)./ew23./ew14)./(n+s);
ayxn = zeros(1,nx*ny);

%ayxw = (exx3.*exx2./ezz2./ew23 - ...
%        exx4.*exx1./ezz1./ew14 + ...
%        e.*(exx4.*exx2 - exx1.*exx3)./ew23./ew14)./(n+s);
ayxs = (ones(1,nx*ny)-exx./ezz)./dx./dy;

%axye = (eyy4.*(1-eyy3./ezz3) - eyy3.*(1-eyy4./ezz4))./ns34./(e+w) - ...
%       2*(eyx1.*eyy2./ezz1.*n.*w./ns21 + ...
%          eyx2.*eyy1./ezz2.*s.*w./ns21 + ...
%          eyx4.*eyy3./ezz4.*n.*e./ns34 + ...
%          eyx3.*eyy4./ezz3.*s.*e./ns34 + ...
%          eyy1.*eyy2.*(1./ezz1-1./ezz2).*w.^2./ns21 + ...
%          eyy3.*eyy4.*(1./ezz4-1./ezz3).*e.*w./ns34)./e./(e+w).^2;
axye = zeros(1,nx*ny);

%  axyw = (eyy2.*(1-eyy1./ezz1) - eyy1.*(1-eyy2./ezz2))./ns21./(e+w) - ...
%         2*(eyx4.*eyy3./ezz4.*n.*e./ns34 + ...
%            eyx3.*eyy4./ezz3.*s.*e./ns34 + ...
%            eyx1.*eyy2./ezz1.*n.*w./ns21 + ...
%            eyx2.*eyy1./ezz2.*s.*w./ns21 + ...
%          eyy4.*eyy3.*(1./ezz3-1./ezz4).*e.^2./ns34 + ...
%          eyy2.*eyy1.*(1./ezz2-1./ezz1).*w.*e./ns21)./w./(e+w).^2;
axyw = (ones(1,nx*ny)-eyy./ezz)./dx./dy;

%ayxn = (exx4.*(1-exx1./ezz1) - exx1.*(1-exx4./ezz4))./ew14./(n+s) - ...
%       2*(exy3.*exx2./ezz3.*e.*s./ew23 + ...
%          exy2.*exx3./ezz2.*w.*s./ew23 + ...
%          exy4.*exx1./ezz4.*e.*n./ew14 + ...
%          exy1.*exx4./ezz1.*w.*n./ew14 + ...
%          exx3.*exx2.*(1./ezz3-1./ezz2).*s.^2./ew23 + ...
%          exx1.*exx4.*(1./ezz4-1./ezz1).*n.*s./ew14)./n./(n+s).^2;
ayxe = (ones(1,nx*ny)-exx./ezzx)./dx./dy;

%ayxs = (exx2.*(1-exx3./ezz3) - exx3.*(1-exx2./ezz2))./ew23./(n+s) - ...
%       2*(exy4.*exx1./ezz4.*e.*n./ew14 + ...
%          exy1.*exx4./ezz1.*w.*n./ew14 + ...
%          exy3.*exx2./ezz3.*e.*s./ew23 + ...
%          exy2.*exx3./ezz2.*w.*s./ew23 + ...
%          exx4.*exx1.*(1./ezz1-1./ezz4).*n.^2./ew14 + ...
%          exx2.*exx3.*(1./ezz2-1./ezz3).*s.*n./ew23)./s./(n+s).^2;
ayxw = zeros(1,nx*ny);

%axyne = +eyy3.*(1-eyy4./ezz4)./(e+w)./ns34;
%axyse = -eyy4.*(1-eyy3./ezz3)./(e+w)./ns34;
%axynw = -eyy2.*(1-eyy1./ezz1)./(e+w)./ns21;
%axysw = +eyy1.*(1-eyy2./ezz2)./(e+w)./ns21;
axyne = zeros(1,nx*ny); %axyse = zeros(1,nx*ny);
axynw = -(ones(1,nx*ny)-eyy./ezzy)./dx./dy; axysw = zeros(1,nx*ny);

%ayxne = +exx1.*(1-exx4./ezz4)./(n+s)./ew14;
%ayxse = -exx2.*(1-exx3./ezz3)./(n+s)./ew23;
%ayxnw = -exx4.*(1-exx1./ezz1)./(n+s)./ew14;
%ayxsw = +exx3.*(1-exx2./ezz2)./(n+s)./ew23;
ayxne = zeros(1,nx*ny); %ayxnw = zeros(1,nx*ny);
ayxse = -(ones(1,nx*ny)-exx./ezzx)./dx./dy; ayxsw = zeros(1,nx*ny);

%axyp = -(axyn + axys + axye + axyw + axyne + axyse + axynw + axysw) ...
%       - k^2.*(w.*(n.*eyx1.*eyy2 + s.*eyx2.*eyy1)./ns21 + ...
%               e.*(s.*eyx3.*eyy4 + n.*eyx4.*eyy3)./ns34)./(e+w);
axyp = -(ones(1,nx*ny)-eyy./ezz)./dx./dy;

%ayxp = -(ayxn + ayxs + ayxe + ayxw + ayxne + ayxse + ayxnw + ayxsw) ...
%       - k^2.*(n.*(w.*exy1.*exx4 + e.*exy4.*exx1)./ew14 + ...
%               s.*(w.*exy2.*exx3 + e.*exy3.*exx2)./ew23)./(n+s);  
ayxp = -(ones(1,nx*ny)-exx./ezz)./dx./dy;

ii = zeros(nx,ny);
ii(:) = (1:nx*ny);

% NORTH boundary

ib = zeros(nx,1);  ib(:) = ii(1:nx,ny);

switch (boundary(1))
  case 'S',   sign = +1;
  case 'A',   sign = -1;
  case '0',   sign = 0;
  otherwise,  
    error('Unrecognized north boundary condition: %s.\n', boundary(1));
end

axxs(ib)  = axxs(ib)  + sign*axxn(ib);
%axxse(ib) = axxse(ib) + sign*axxne(ib);
%axxsw(ib) = axxsw(ib) + sign*axxnw(ib);
ayxs(ib)  = ayxs(ib)  + sign*ayxn(ib);
%ayxse(ib) = ayxse(ib) + sign*ayxne(ib);
%ayxsw(ib) = ayxsw(ib) + sign*ayxnw(ib);
ayys(ib)  = ayys(ib)  - sign*ayyn(ib);
%ayyse(ib) = ayyse(ib) - sign*ayyne(ib);
%ayysw(ib) = ayysw(ib) - sign*ayynw(ib);
axys(ib)  = axys(ib)  - sign*axyn(ib);
%axyse(ib) = axyse(ib) - sign*axyne(ib);
axysw(ib) = axysw(ib) - sign*axynw(ib);

% SOUTH boundary

ib = zeros(nx,1);  ib(:) = ii(1:nx,1);

switch (boundary(2))
  case 'S',   sign = +1;
  case 'A',   sign = -1;
  case '0',   sign = 0;
  otherwise,  
    error('Unrecognized south boundary condition: %s.\n', boundary(2));
end

axxn(ib)  = axxn(ib)  + sign*axxs(ib);
%axxne(ib) = axxne(ib) + sign*axxse(ib);
%axxnw(ib) = axxnw(ib) + sign*axxsw(ib);
ayxn(ib)  = ayxn(ib)  + sign*ayxs(ib);
ayxne(ib) = ayxne(ib) + sign*ayxse(ib);
%ayxnw(ib) = ayxnw(ib) + sign*ayxsw(ib);
ayyn(ib)  = ayyn(ib)  - sign*ayys(ib);
%ayyne(ib) = ayyne(ib) - sign*ayyse(ib);
%ayynw(ib) = ayynw(ib) - sign*ayysw(ib);
axyn(ib)  = axyn(ib)  - sign*axys(ib);
%axyne(ib) = axyne(ib) - sign*axyse(ib);
%axynw(ib) = axynw(ib) - sign*axysw(ib);

% EAST boundary

ib = zeros(1,ny);  ib(:) = ii(nx,1:ny);

switch (boundary(3))
  case 'S',   sign = +1;
  case 'A',   sign = -1;
  case '0',   sign = 0;
  otherwise,  
    error('Unrecognized east boundary condition: %s.\n', boundary(3));
end

axxw(ib)  = axxw(ib)  + sign*axxe(ib);
%axxnw(ib) = axxnw(ib) + sign*axxne(ib);
%axxsw(ib) = axxsw(ib) + sign*axxse(ib);
ayxw(ib)  = ayxw(ib)  + sign*ayxe(ib);
%ayxnw(ib) = ayxnw(ib) + sign*ayxne(ib);
ayxsw(ib) = ayxsw(ib) + sign*ayxse(ib);
ayyw(ib)  = ayyw(ib)  - sign*ayye(ib);
%ayynw(ib) = ayynw(ib) - sign*ayyne(ib);
%ayysw(ib) = ayysw(ib) - sign*ayyse(ib);
axyw(ib)  = axyw(ib)  - sign*axye(ib);
%axynw(ib) = axynw(ib) - sign*axyne(ib);
%axysw(ib) = axysw(ib) - sign*axyse(ib);

% WEST boundary

ib = zeros(1,ny);  ib(:) = ii(1,1:ny);

switch (boundary(4))
  case 'S',   sign = +1;
  case 'A',   sign = -1;
  case '0',   sign = 0;
  otherwise,  
    error('Unrecognized west boundary condition: %s.\n', boundary(4));
end

axxe(ib)  = axxe(ib)  + sign*axxw(ib);
%axxne(ib) = axxne(ib) + sign*axxnw(ib);
%axxse(ib) = axxse(ib) + sign*axxsw(ib);
ayxe(ib)  = ayxe(ib)  + sign*ayxw(ib);
%ayxne(ib) = ayxne(ib) + sign*ayxnw(ib);
%ayxse(ib) = ayxse(ib) + sign*ayxsw(ib);
ayye(ib)  = ayye(ib)  - sign*ayyw(ib);
%ayyne(ib) = ayyne(ib) - sign*ayynw(ib);
%ayyse(ib) = ayyse(ib) - sign*ayysw(ib);
axye(ib)  = axye(ib)  - sign*axyw(ib);
axyne(ib) = axyne(ib) - sign*axynw(ib);
%axyse(ib) = axyse(ib) - sign*axysw(ib);

% Assemble sparse matrix

iall = zeros(1,nx*ny);          iall(:) = ii;
is = zeros(1,nx*(ny-1));        is(:) = ii(1:nx,1:(ny-1));
in = zeros(1,nx*(ny-1));        in(:) = ii(1:nx,2:ny);
ie = zeros(1,(nx-1)*ny);        ie(:) = ii(2:nx,1:ny);
iw = zeros(1,(nx-1)*ny);        iw(:) = ii(1:(nx-1),1:ny);
ine = zeros(1,(nx-1)*(ny-1));   ine(:) = ii(2:nx, 2:ny);
ise = zeros(1,(nx-1)*(ny-1));   ise(:) = ii(2:nx, 1:(ny-1));
isw = zeros(1,(nx-1)*(ny-1));   isw(:) = ii(1:(nx-1), 1:(ny-1));
inw = zeros(1,(nx-1)*(ny-1));   inw(:) = ii(1:(nx-1), 2:ny);

Axx = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[axxp(iall),axxe(iw),axxw(ie),axxn(is),axxs(in)]);

Axy = sparse ([iall,iw,ie,is,in,ine,ise,isw], ...
	[iall,ie,iw,in,is,isw,inw,ine], ...
	[axyp(iall),axye(iw),axyw(ie),axyn(is),axys(in), ...
     axysw(ine),axynw(ise),axyne(isw)]);

Ayx = sparse ([iall,iw,ie,is,in,ine,isw,inw], ...
	[iall,ie,iw,in,is,isw,ine,ise], ...
	[ayxp(iall),ayxe(iw),ayxw(ie),ayxn(is),ayxs(in), ...
     ayxsw(ine),ayxne(isw),ayxse(inw)]);

Ayy = sparse ([iall,iw,ie,is,in], ...
	[iall,ie,iw,in,is], ...
	[ayyp(iall),ayye(iw),ayyw(ie),ayyn(is),ayys(in)]);

A = [[Axx Axy];[Ayx Ayy]];

fprintf(1,'nnz(A) = %d\n',nnz(A));
fprintf(1,'density(A) = %d\n',nnz(A)/prod(size(A)));

shift = (guess*k)^2;
options.tol = 1e-8;
options.disp = 0;						% suppress output

clear Axx Axy Ayx Ayy ...
    axxnw axxne axxne ...
    axxw  axxp  axxe ...
    axxsw axxse axxse ...
    axynw axyne axyne ...
    axyw  axyp  axye ...
    axysw axyse axyse ...
    ayynw ayyne ayyne ...
    ayyw  ayyp  ayye ...
    ayysw ayyse ayyse ...
    ayxnw ayxne ayxne ...
    ayxw  ayxp  ayxe ...
    ayxsw ayxse ayxse ...
    iall ie iw in iw ...
    isw inw ine ise ...
    exx1 exx2 exx3 exx4 ...
    exy1 exy2 exy3 exy4 ...
    eyx1 eyx2 eyx3 eyx4 ...
    eyy1 eyy2 eyy3 eyy4 ...
    ezz1 ezz2 ezz3 ezz4 ...
    ns21 ns34 ew14 ew23;

[v,d] = eigs(A,speye(size(A)),nmodes,shift,options);
neff = lambda*sqrt(diag(d))/(2*pi);
clear d;

phix = zeros(nx,ny,nmodes);
phiy = zeros(nx,ny,nmodes);
%temp = zeros(nx,2*ny);

% Normalize modes

temp = zeros(nx*ny,2);

for kk = 1:nmodes;
  temp(:) = v(:,kk);
  [mag,ii] = max(sqrt(sum(abs(temp).^2,2)));
  if abs(temp(ii,1)) > abs(temp(ii,2)),
    jj = 1;
  else 
    jj = 2;
  end
  mag = mag*temp(ii,jj)/abs(temp(ii,jj));
  temp = temp/mag;
  phix(:,:,kk) = reshape(temp(:,1),nx,ny);
  phiy(:,:,kk) = reshape(temp(:,2),nx,ny);

end;

phiz = zeros(nx-1,ny-1,nmodes);
psix = zeros(nx,ny-1,nmodes);
psiy = zeros(nx-1,ny,nmodes);
psiz = zeros(nx-1,ny-1,nmodes);

for kk = 1:nmodes;
  phiz(:,:,kk) = -1.0./(neff(kk).*k).*((phix(2:nx,1:ny-1,kk)-phix(1:nx-1,1:ny-1,kk))./dx + ...
                  (phiy(1:nx-1,2:ny,kk)-phiy(1:nx-1,1:ny-1,kk))./dy);
  phiz = [phiz(:,:,kk),phiz(:,ny-1,kk)];
  phiz = [phiz(:,:,kk);phiz(nx-1,:,kk)];

  psix(:,:,kk) = 1.0./(k.*epsxx(1:nx,2:ny)).*((phiz(1:nx,2:ny,kk)-phiz(1:nx,1:ny-1,kk))./dy);
  psix = [psix(:,1,kk),psix(:,:,kk)];%,psix(:,ny-1,kk)];
  psix(:,:,kk) = psix(:,:,kk) + neff(kk)./epsxx(:,:).*phiy(:,:,kk);

  psiy(:,:,kk) = -1.0./(k.*epsyy(2:nx,1:ny)).*((phiz(2:nx,1:ny,kk)-phiz(1:nx-1,1:ny,kk))./dx);
  psiy = [psiy(1,:,kk);psiy(:,:,kk)];%;psiy(nx-1,:,kk)];
  psiy(:,:,kk) = psiy(:,:,kk) - neff(kk)./epsyy(:,:).*phix(:,:,kk);

  psiz(:,:,kk) = -1.0./(k.*epszz(2:nx,2:ny)).*((phiy(2:nx,2:ny,kk) - phiy(1:nx-1,2:ny,kk))./dx - ...
                  (phix(2:nx,2:ny,kk) - phix(2:nx,1:ny-1,kk))./dy);
  psiz = [psiz(:,1,kk),psiz(:,:,kk)];%,psiz(:,ny-1,kk)];
  psiz = [psiz(1,:,kk);psiz(:,:,kk)];%;psiz(nx-1,:,kk)];
end;
fprintf (1,'done (cputime = %7.3f)\n\n', cputime-t);

mag = max(max(psix));
phix = phix./mag; phiy = phiy./mag; phiz = phiz./mag;
psix = psix./mag; psiy = psiy./mag; psiz = psiz./mag;

return;

