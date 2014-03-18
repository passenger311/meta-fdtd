function output = SL_NRI(omega)

w_p     =   2*pi*10.0013045e9;
gamma_d =   4.45146831e8;
n_d     =   1-((w_p^2)/((omega^2)+complex(0,(omega*gamma_d))));
eps_d   =   n_d;
mu_d    =   eps_d;

output  =   [eps_d,mu_d];