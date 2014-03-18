% Complex frequency drude derivative


function out = cf_drude_K(w)

eps_inf     =   5;
wp          =   1.4433e16;
gamma       =   0;%6.0771e13;

deps        =   (2*w+(complex(0,gamma)))*((wp^2)/((((w^2)+complex(0,1)*gamma*w))^2));
eps         =   eps_inf-((wp^2)/((w^2)+complex(0,1)*gamma*w));

out     =   [eps deps];