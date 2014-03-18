
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

function par = cgd(p_par)



% Calculate other photon parameters

lambda  =   p_par(1);               % Free space wavelength (m)
omega   =   p_par(2);               % Free space angular frequency (rad/s)
freq    =   p_par(3);               % Free space frequency (1/s)
k0      =   p_par(4);               % Free space wavevector (1/m)
E       =   p_par(5);               % Free space photon energy (eV)

% Step 1: Define cover and substrate layers

epsc    =   4.26^2;                 % Permittivity of cover layer
epss    =   complex(0.583,10.1)^2;  % Permittivity of substrate layer

muc     =   1;                      % Permeability of cover layer
mus     =   1;                      % Permeability of substrate layer

% Step 2: Define internal layers (1:N)

eps_layers  =   [1.44^2];           % Permittivity of layers (1:N)
mu_layers   =   [1];                % Permeability of layers (1:N)
d_layers    =   [2e-9];             % Thickness of layers (1:N) (m)


% Build structure parameter arrays

eps_struc   =   [epsc eps_layers epss];     % Structure permittivity array
mu_struc    =   [muc mu_layers mus];        % Structure permeabillity array
d_struc     =   [0 0 d_layers];             % Structure thickness array (zeros included to match length for output)     
n_struc     =   sqrt(eps_struc.*mu_struc);  % Structure permittivity array
z_struc     =   sqrt(mu_struc./eps_struc);  % Structure impedance array


% Build output matrix 

par(1,:)  =   eps_struc;
par(2,:)  =   mu_struc;
par(3,:)  =   n_struc;
par(4,:)  =   z_struc;
par(5,:)  =   d_struc;