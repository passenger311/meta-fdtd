function output = E_gain(omega)

eps_inf     =   1.001;
delta_eps   =   0.0053;
wl          =   2*pi*370e12;
gammal      =   0.1e14;%1e14;

eps_l       =   eps_inf+((delta_eps*(wl^2))/((wl^2)+complex(0,(2*gammal*omega))-(omega^2)));

output      =   [eps_l,1];