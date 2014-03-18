function output = E_NRI(omega)

w_p     =   2*pi*893.8e12;
gamma_d =   0.27e12;
n_d     =   1-((w_p^2)/((omega^2)+complex(0,(omega*gamma_d))));
eps_d   =   n_d;
mu_d    =   eps_d;

output  =   [eps_d,mu_d];