function sol = dD_dk(w,b,pola,c,name,d_val)

fh          =   str2func(name);
inputs      =   fh(w,b,pola,c);
kappa_val   =   inputs{1};
m_val       =   inputs{3};

%d_val       =   inputs{5};
gamma1      =   inputs{6};

gammaN      =   inputs{8};



% Initialise phi and psi
phi1    =   1;
phi2    =   -gamma1/m_val(1);
phi     =   [phi1;phi2];

psi1    =   0;
psi2    =   -1/(2*m_val(1));
psi     =   [psi1;psi2];

N       =   length(m_val);

for layer = 2:N-1

    d       =   d_val(layer-1);
    m       =   m_val(layer);
    kappa   =   kappa_val(layer-1);

    M       =   [cos(kappa*d) (-(m/kappa)*sin(kappa*d))];
    M       =   [M; (kappa/m)*sin(kappa*d) (cos(kappa*d))];

    DM0      =  (d/(2*kappa))*sin(kappa*d);
    DM1      =  ((-m/(2*(kappa^3)))*sin(kappa*d))+(m*d/(2*(kappa^2)))*cos(kappa*d);
    DM2      =  ((-1/(2*kappa*m))*sin(kappa*d))-((d/(2*m))*cos(kappa*d));
    DM3     =   (d/(2*kappa))*sin(kappa*d);

    DM      =   [DM0 DM1; DM2 DM3];

    psi     =   M*psi + gamma1*DM*phi;
    phi     =   M*phi;

end

H   =   phi(2)-(gammaN/m_val(end))*phi(1);

DH  =   2*b*( (psi(2)/gamma1) - (phi(1)/(2*gammaN*m_val(end))) - (gammaN/(gamma1*m_val(end)))*psi(1));

sol = [H DH];



