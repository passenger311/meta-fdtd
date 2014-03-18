function vg = g_vel(z_in,s_par,k0,name,pola,R,S,meth,nrs,d_struc)

p_par   =   photon_par([],[],[],k0,[]);
c_lambda =   p_par(1);
h       =   c_lambda/1e10;
w       =   exp(2*z_in)+((R^2)*exp(-2*z_in))+S;
n_c     =   sqrt(w);
low     =   c_lambda-h;

par1    =   var(low,1,1,1,[],'lambda',1,s_par,p_par,0,name,1);

% Update other parameters
k0      =   par1(1,2);
s_par   =   par1(3:end,:);

eps_struc   =   s_par(1,:);
mu_struc    =   s_par(2,:);
n_struc     =   s_par(3,:);


nc          =   n_struc(1);
ns          =   n_struc(end);
R           =   ((ns^2) - (nc^2))/4;       % Coefficient used in mapping
S           =   ((ns^2) + (nc^2))/2;        % Coefficient used in mapping

if pola == 1
    m_struc = s_par(2,:);
elseif pola == 0;
    m_struc = s_par(1,:);
end

z_low       =   NR(eps_struc,d_struc,m_struc,mu_struc,R,S,z_in,k0,nrs);
w           =   exp(2*z_low)+((R^2)*exp(-2*z_low))+S;
n_l         =   sqrt(w);
high        =   c_lambda+h;


par1    =   var(high,1,1,1,[],'lambda',1,s_par,p_par,0,name,1);

% Update other parameters
k0      =   par1(1,2);
s_par   =   par1(3:end,:);

eps_struc   =   s_par(1,:);
mu_struc    =   s_par(2,:);
n_struc     =   s_par(3,:);


nc          =   n_struc(1);
ns          =   n_struc(end);
R           =   ((ns^2) - (nc^2))/4;        % Coefficient used in mapping
S           =   ((ns^2) + (nc^2))/2;        % Coefficient used in mapping

if pola == 1
    m_struc = s_par(2,:);
elseif pola == 0;
    m_struc = s_par(1,:);
end

z_high      =   NR(eps_struc,d_struc,m_struc,mu_struc,R,S,z_in,k0,nrs);
w           =   exp(2*z_high)+((R^2)*exp(-2*z_high))+S;
n_h         =   sqrt(w);

if meth == 1
    dn_dl   =   (real(n_h)-real(n_l))/(2*h);
elseif meth == 2
    dn_dl   =   (real(n_h)-real(n_c))/h;
else
    dn_dl   =   (real(n_c)-real(n_l))/h;
end

vg      =   1/(real(n_c)-(dn_dl*c_lambda));
