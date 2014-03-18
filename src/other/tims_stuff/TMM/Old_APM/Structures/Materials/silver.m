function output = silver(w)

eps_inf =   1.17152;
wd      =   1.39604e16;
gammad  =   12.6126;
dl1     =   2.23994;
wl1     =   8.25718e15;
gammal1 =   1.95614e14;
dl2     =   0.222651;
wl2     =   3.05707e15;
gammal2 =   8.52675e14;

i       =   complex(0,1);
eps_s   =   eps_inf - (wd^2)/((w^2)-i*gammad*w);
eps_s   =   eps_s - (dl1*(wl1^2))/((w^2)-(i*2*gammal1*w)-(wl1^2));
eps_s   =   eps_s - (dl2*(wl2^2))/((w^2)-(i*2*gammal2*w)-(wl2^2));

output  =   eps_s;