function output = ITO(w)

wp          =   3.13e15;
gamma       =   1.07e14;
eps_inf     =   4;
eps         =   eps_inf-((wp^2)./((w.^2)+1i*gamma*w));


output = (eps);