function output = drude(eps_inf,wp,gamma,w)

eps         =   eps_inf-((wp^2)/((w^2)+1i*gamma*w));

output = eps;