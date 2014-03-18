function output = lossy_drude(w)

wp          =   1;
gamma       =   0;
eps_inf     =   1;
eps         =   eps_inf-((wp^2)/((w^2)+complex(0,gamma*w)));


output = eps;