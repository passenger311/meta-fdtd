
% -------------------------------------------------------------------------
% |                                                                       |
% |                            Function Struc                             |
% |                                                                       |
% | In this function you define the layered structure you want to solve.  |
% | The structure is always of the form shown below.  When solving for    |
% | reflection and transmission coefficients light is always incident     |
% | from the cover layer.                                                 |
% |                                                                       |               
% |                                                                       |
% |                      |                        |                       |                         
% |          Cover layer | .... Layers (1:N) .... | Substrate layer       |
% |             (epsc)   |          eps(1:N)      |     (epss)            |
% |             (muc)    |           mu(1:N)      |     (mus)             |
% |                      |            d(1:N)      |                       |
% |                                                                       |
% -------------------------------------------------------------------------

function par = Test1(p_par)



% Calculate other photon parameters

lambda  =   p_par(1);               % Free space wavelength (m)
omega   =   p_par(2);               % Free space angular frequency (rad/s)
freq    =   p_par(3);               % Free space frequency (1/s)
k0      =   p_par(4);               % Free space wavevector (1/m)
E       =   p_par(5);               % Free space photon energy (eV)

% Step 1: Define cover and substrate layers


ag      =   (-143.49-1i*9.52);
siO2    =   (2.1025);
au      =   (-95.92-1i*10.97);
eps0    =   1.45^2;
eps1    =   1.39^2;




d_layers    =   [1.57 0.78 1.57].*lambda;             % Thickness of layers (1:N) (m)


% Build structure parameter arrays

eps_struc   =   [eps0 eps1 eps0 eps1 eps0];                % Structure permittivity array
mu_struc    =   [1 1 1 1 1];                      % Structure permeabillity array
d_struc     =   [0 0 d_layers];             % Structure thickness array (zeros included to match length for output)     
n_struc     =   sqrt(eps_struc.*mu_struc);  % Structure permittivity array
z_struc     =   sqrt(mu_struc./eps_struc);  % Structure impedance array


% Build output matrix 

par(1,:)  =   eps_struc;
par(2,:)  =   mu_struc;
par(3,:)  =   n_struc;
par(4,:)  =   z_struc;
par(5,:)  =   d_struc;