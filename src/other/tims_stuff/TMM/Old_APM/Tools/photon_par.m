
% ------------------------------------------------------------------------
% |                                                                      |
% |                          Function photon_par                         |
% |                                                                      |
% | This function takes one of 5 photon parameters (Wavelength, Angular  |
% | frequency, Frequency, Angular wavenumber and photon energy) as an    |
% | input then calculates and outputs an array (par) containing all 5    |
% | parameters. This is mainly for simplifying the input stage.          |
% |                                                                      |
% ------------------------------------------------------------------------


function par = photon_par(lambda,omega,freq,k0,E)

% Constants 
c   =   299792458;              % Speed of light (m/s)
h   =   4.13566733e-15;         % Plancks constant (eV/s)

% Check inputs to find initial parameter to use to calculate others
check(1)    =   isempty(lambda);
check(2)    =   isempty(omega);
check(3)    =   isempty(freq);
check(4)    =   isempty(k0);
check(5)    =   isempty(E);

ind         =   find(check == 0);

% If more than one input parameter is defined print error message and quit
if isempty(ind) == 1 || length(ind) > 1
    fprintf('Please only define one photon parameter\n\n')
    return
end

% From initial parameter calculate other photon parameters
if ind == 1                         % Input is wavelength (lambda)
    omega   =   (2*pi*c)/lambda;
    freq    =   c/lambda;
    k0      =   (2*pi)/lambda;
    E       =   (h*c)/lambda;
elseif ind == 2                     % Input is angular frequency (omega)
    lambda  =   (2*pi*c)/omega;
    freq    =   c/lambda;
    k0      =   (2*pi)/lambda;
    E       =   (h*c)/lambda;
elseif ind == 3                     % Input is frequency (freq)
    lambda  =   c/freq;
    omega   =   (2*pi*c)/lambda;
    k0      =   (2*pi)/lambda;
    E       =   (h*c)/lambda;
elseif ind == 4                     % Input is angular wavenumber (k0)
    lambda  =   (2*pi)/k0;
    freq    =   c/lambda;
    omega   =   (2*pi*c)/lambda;
    E       =   (h*c)/lambda;
elseif ind == 5                     % Input is photon energy (E)
    lambda  =   (h*c)/E;
    omega   =   (2*pi*c)/lambda;
    freq    =   c/lambda;
    k0      =   (2*pi)/lambda;
end

% Bundle photon parameters into 'par' to output
par(1)  =   lambda;
par(2)  =   omega;
par(3)  =   freq;
par(4)  =   k0;
par(5)  =   E;
    



