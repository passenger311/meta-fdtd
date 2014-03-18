% Complex frequency drude derivative

function out = cf_drude(w)

eps_inf     =   1;
wp          =   1;
gamma       =   0;%1e-2;

deps        =   (2*w+(complex(0,gamma)))*((wp^2)/((((w^2)+complex(0,1)*gamma*w))^2));
eps         =   eps_inf-((wp^2)/((w^2)+complex(0,1)*gamma*w));

out     =   [eps deps];