function in_par = input_parameters


% -------------------------- General parameters --------------------------

% Structure parameters
name    =   'zeph';                 % Name of structure

% Photon parameters

lambda  =   [1.55e-6];              % Free space wavelength (m)
omega   =   [];                     % Free space angular frequency (rad/s)
freq    =   [];                   % Free space frequency (1/s)
k0      =   [];                     % Free space wavevector (rad/m)
E       =   [];                     % Free space photon energy (eV)

% -------------------------- Testing parameters -------------------------- 

% Tests to run
m_in    =   'y';                  % Perform initial scan for modes (Yes/No)
rt_in   =   'n';                   % Perform initial reflection calculation (Yes/No)
m_scan  =   1;                      % Track modes (number of runs)
rt_scan =   0;                      % Track reflection trans (number of runs)

% Save parameters
prec    =   16;                     % Output significant figures
save_i  =   'Yes';                  % Save initial results to file (Yes/No)

% Modal analysis specific parameters

% Contour parameters
re0     =  -4;                     % Min real val of complex plane
re1     =   5;                      % Max real val of complex plane
im0     =   -1.1*pi;                % Min imag val of complex plane
im1     =   1.2*pi;                 % Max imag val of complex plane

% Other parameters
pola    =   'TM';                   % Mode polarisation (TE or TM)
nrs     =   100;                     % Number of Newton Raphson steps
vg      =   'n';                   % Calculate group velocity?
meth    =   1;                      % Method

% Reflection transmission specific parameters

theta_i =   0;                     % Initial angle of incidence (degrees)


% --------------------- Bundle parameters and output ---------------------

% Calculate other photon parameters and bundle
p_par   =   photon_par(lambda,omega,freq,k0,E);

% Bundle contour parameters
cont    =   [re0 re1 im0 im1];

% Digitise other parameters
m_in    =   strcmpi(m_in(1),'y');   % Digitise m_in
rt_in   =   strcmpi(rt_in(1),'y');  % Digitise rt_in
save_i  =   strcmpi(save_i(1),'y'); % Digitise save choice
pola    =   strcmpi(pola(end),'e'); % Digitise polarisation choice
vg      =   strcmpi(vg(1),'y');     % Digitise group velocity choice

% Bundle other parameters
other   =   [m_in rt_in m_scan rt_scan prec save_i]; 

% Create output array
in_par = [p_par other cont pola nrs theta_i vg meth];
in_par = {name in_par};
