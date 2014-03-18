% Complex frequency chl

function out = cf_chl(w,b,pola,c)

eps         =   cf_drude(w);
epsp        =   eps(1);
depsp       =   eps(2);
eps_struc   =   [1.44^2 4.27^2 epsp];
mu_struc    =   [1 1 1];
n2_struc    =   eps_struc.*mu_struc;
n_struc     =   sqrt(n2_struc);
lambdap     =   (2*pi*c);
d_struc     =   [0.02*lambdap];

deps_struc  =   [0 0 depsp];
dmu_struc   =   [0 0 0];

gamma1      =   sqrt((b^2)-((w/c)^2)*n2_struc(1));

n1          =   n_struc(1);
dn1         =   deps_struc(1);
dgamma1     =   -w*n1*(w*dn1+n1)/(gamma1*(c^2));
gammaN      =   sqrt((b^2)-((w/c)^2)*epsp);
nN          =   n_struc(end);
dnN         =   deps_struc(end)/(2*nN);

dgammaN     =   -w*nN*(w*dnN+nN)/(gammaN*(c^2));

kappa_struc =   sqrt((((w/c)^2)*n2_struc(2))-(b^2));
dkappa_struc    =   w*(n2_struc(2))/(kappa_struc*(c^2));

if pola == 1
    m_struc     =   mu_struc;
    dm_struc    =   dmu_struc;
else
    m_struc     =   eps_struc;
    dm_struc    =   deps_struc;
end
    


out = {kappa_struc,dkappa_struc,m_struc,dm_struc,d_struc,gamma1,dgamma1,gammaN,dgammaN};



