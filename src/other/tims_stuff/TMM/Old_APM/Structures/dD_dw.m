function sol = dD_dw(w,b,pola,c,name,d_val)

fh          =   str2func(name);
inputs      =   fh(w,b,pola,c);
kappa_val   =   inputs{1};
dkappa_val  =   inputs{2};
m_val       =   inputs{3};
dm_val      =   inputs{4};
%d_val       =   inputs{5};
gamma1      =   inputs{6};
dgamma1     =   inputs{7};
gammaN      =   inputs{8};
dgammaN     =   inputs{9};


% Initialise phi and psi
phi1    =   1;
phi2    =   -gamma1/m_val(1);
phi     =   [phi1;phi2];

psi1    =   0;
psi2    =   ((gamma1*dm_val(1))-(dgamma1*m_val(1)))/(m_val(1)^2);
psi     =   [psi1;psi2];

N       =   length(m_val);

for layer = 2:N-1

    d       =   d_val(layer-1);
    m       =   m_val(layer);
    dm      =   dm_val(layer);
    kappa   =   kappa_val(layer-1);
    d_k     =   dkappa_val(layer-1);

    M       =   [cos(kappa*d) (-(m/kappa)*sin(kappa*d))];
    M       =   [M; (kappa/m)*sin(kappa*d) (cos(kappa*d))];

    DM0     =   -d*d_k*sin(d*kappa);
    DM1     =   (((m*d_k)-(dm*kappa))/(kappa^2))*sin(kappa*d);
    DM1     =   DM1 - ((m*d*d_k/kappa)*cos(kappa*d));
    DM2     =   (((m*d_k)-(dm*kappa))/(m^2))*sin(kappa*d);
    DM2     =   DM2 + ((d*kappa*d_k/m)*cos(kappa*d));  

    DM      =   [DM0 DM1; DM2 DM0];

    psi     =   M*psi + DM*phi;
    phi     =   M*phi;

end

H   =   phi(2)-(gammaN/m_val(end))*phi(1);

DH  =   (psi(2)) + (((gammaN*dm_val(end))-(dgammaN*m_val(end)))/(m_val(end)^2))*phi(1)-((gammaN*psi(1))/(m_val(end)));

sol = [H DH];



