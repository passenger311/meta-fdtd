
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

function par = J_3_layer(p_par)



% Calculate other photon parameters

lambda  =   p_par(1);               % Free space wavelength (m)
omega   =   p_par(2);               % Free space angular frequency (rad/s)
freq    =   p_par(3);               % Free space frequency (1/s)
k0      =   p_par(4);               % Free space wavevector (1/m)
E       =   p_par(5);               % Free space photon energy (eV)

% Step 1: Define cover and substrate layers
c           =   299792458;
air         =   1;
Si          =   3.47^2;
SiO2        =   1.44^2;

d_layers    =   [80 200]*1e-9;             % Thickness of layers (1:N) (m)


% Build structure parameter arrays

eps_struc   =   [air SiO2 Si SiO2];                % Structure permittivity array
mu_struc    =   [1 1 1 1];                      % Structure permeabillity array
d_struc     =   [0 0 d_layers];             % Structure thickness array (zeros included to match length for output)     
n_struc     =   sqrt(eps_struc.*mu_struc);  % Structure permittivity array
z_struc     =   sqrt(mu_struc./eps_struc);  % Structure impedance array


% Build output matrix 

par(1,:)  =   eps_struc;
par(2,:)  =   mu_struc;
par(3,:)  =   n_struc;
par(4,:)  =   z_struc;
par(5,:)  =   d_struc;