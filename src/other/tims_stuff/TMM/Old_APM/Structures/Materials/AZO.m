function output = AZO(w)

wp          =   2.9e15;
gamma       =   1.6e14;
eps_inf     =   4;
eps         =   eps_inf-((wp^2)/((w^2)+complex(0,gamma*w)));


output = real(eps);