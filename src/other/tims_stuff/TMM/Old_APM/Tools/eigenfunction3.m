function sol = eigenfunction3(eps_val,d_val,m_val,mu_val,b,k0)


n2_val      =   eps_val.*mu_val;
N           =   length(eps_val);
gamma1  =   sqrt(b^2-(k0^2)*n2_val(1));
gammaN  =   sqrt(b^2-(k0^2)*n2_val(end));


phi         =   [1; -gamma1/m_val(1)];
psi         =   [0;(-1/(2*m_val(1)))];

for layer = 2:N-1

    d       =   d_val(layer-1);
    m       =   m_val(layer);
    n2      =   n2_val(layer);
    kappa   =   sqrt((k0^2)*n2-b^2);

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



