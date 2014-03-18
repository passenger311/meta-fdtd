function output = lossy_drude_K(w)

eps_inf     =   5;
wp          =   1.4433e16;
gamma       =   6.0771e13;


eps =   eps_inf-((wp^2)/((w^2)+complex(0,1)*gamma*w));



output = eps;