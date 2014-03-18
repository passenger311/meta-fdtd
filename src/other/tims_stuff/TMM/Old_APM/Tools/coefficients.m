function c = coefficients(n_struc,z_struc,d_struc,theta_i,k0,~)

% n_struc is the refractive indicies of the structure starting from the
% incident side.
% z_struc is the impedances of the layers.
% d_struc is the thicknesses of the layers in meters.
% theta_i is the incident angle from the normal.
% k0 is the incident freespace wavenumber.

c1      =   [0;1];
c_par   =   c1;
c_per   =   c1;
nsinc   =   n_struc(1)*sin(theta_i*pi/180);
%nsinc   =   11.1686e6/k0;

% Flip structure around

n_struc = fliplr(n_struc);
z_struc = fliplr(z_struc);
d_struc = fliplr(d_struc);
d_struc = [d_struc 0];


% Calculate amplitude coefficients

for loop2 = 1:length(n_struc)-1
    
    n1      =   n_struc(loop2);
    n2      =   n_struc(loop2+1);
    stheta1 =   nsinc/n1;
    stheta2 =   nsinc/n2;
    ctheta1 =   sqrt(1-stheta1^2);
    ctheta2 =   sqrt(1-stheta2^2);
    z1      =   z_struc(loop2);
    z2      =   z_struc(loop2+1);
    d2      =   d_struc(loop2);
    
    a       =   ctheta1/ctheta2;
    b       =   z2/z1;
    s       =   n2*k0*ctheta2;
    e_pos   =   exp(complex(0,-1)*s*d2);
    e_neg   =   exp(complex(0,1)*s*d2);
    
    M_par0  =   ((a+b)*e_neg);
    M_par1  =   ((a-b)*e_neg);
    M_par2  =   ((a-b)*e_pos);
    M_par3  =   ((a+b)*e_pos);
    
    M_par   =   0.5.*[M_par0 M_par1;M_par2 M_par3];
    
    M_per0  =   (1+(a*b))*e_neg;
    M_per1  =   (1-(a*b))*e_neg;
    M_per2  =   (1-(a*b))*e_pos;
    M_per3  =   (1+(a*b))*e_pos;
    
    M_per   =   0.5.*[M_per0 M_per1;M_per2 M_per3];
    
    c_par   =   M_par*c_par;
    c_per   =   M_per*c_per;
    
end

% Calculate Power coefficients

sthetac =   nsinc/n_struc(1);
sthetas =   nsinc/n_struc(end);

cthetac =   sqrt(1-sthetac^2);
cthetas =   sqrt(1-sthetas^2);

zs      =   z_struc(end);
zc      =   z_struc(1);
inc_par =   (c_par(2));
tau_par =   1/inc_par;
ref_par =   c_par(1)/inc_par;
R_par   =   abs(ref_par)^2;


T_par   =   (zs/zc)*(cthetac/cthetas)*abs(tau_par)^2;

inc_per =   c_per(2);
tau_per =   1/inc_per;
ref_per =   c_per(1)/inc_per;
R_per   =   abs(ref_per)^2;


T_per   =   (zs/zc)*(cthetac/cthetas)*abs(tau_per)^2;

%c = [inc_par R_par;inc_per R_per];
 c = [T_par R_par;T_per R_per];