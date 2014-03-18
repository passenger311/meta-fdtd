function sol = eigenfunction(eps_val,d_val,m_val,mu_val,R,S,z,k0)

gammac  =   k0*(exp(z)+R*exp(-z));
gammas  =   k0*(exp(z)-R*exp(-z));
w       =   exp(2*z)+((R^2)*exp(-2*z))+S;

phi         =   [1; -gammac/m_val(1)];
psi         =   [0;(-1/(2*m_val(1)))];

n2_val      =   eps_val.*mu_val;
N           =   length(eps_val);

for layer = 2:N-1

    d       =   d_val(layer-1);
    m       =   m_val(layer);
    n2      =   n2_val(layer);
    kappa   =   k0*sqrt(n2-w);

    M       =   [cos(kappa*d) (-(m/kappa)*sin(kappa*d))];
    M       =   [M; (kappa/m)*sin(kappa*d) (cos(kappa*d))];

    DM0      =  (d/(2*kappa))*sin(kappa*d);
    DM1      =  ((-m/(2*(kappa^3)))*sin(kappa*d))+(m*d/(2*(kappa^2)))*cos(kappa*d);
    DM2      =  ((-1/(2*kappa*m))*sin(kappa*d))-((d/(2*m))*cos(kappa*d));
    DM3     =   (d/(2*kappa))*sin(kappa*d);

    DM      =   [DM0 DM1; DM2 DM3];

    psi     =   M*psi + gammac*DM*phi;
    phi     =   M*phi;

end

H   =   phi(2)-(gammas/m_val(end))*phi(1);

DH  =   2*(gammas*psi(2) - (gammac/(2*m_val(end)))*phi(1) - ((gammas^2)/m_val(end))*psi(1));

sol = [H DH];