
% -----------------------------------------------------------------------
% |                                                                     |
% |                  1D planar multilayer mode solver                   |
% |                                                                     |
% |  This program formulates the eigenvalue equation for a multilayer   |
% |  waveguide using the transfer matrix method (TMM).  Conformal       |
% |  mapping is used to unwrap the four-sheeted Reimann surface of the  |
% |  eigenfunction.  The zeros of the resulting function are found      |
% |  using the argument principle method (APM).  A selected mode is     |
% |  then tracked using the Newton Raphson method while parameters      |
% |  such as wavelength or layer thickness are varied.                  |
% |                                                                     |
% -----------------------------------------------------------------------


close all;
clear all;
%clc
a       =   cd;
ind1    =   strfind(a, 'Final_program');
a       =   a(1:ind1+12);
addpath(genpath(a))
addpath 'Inputs'
addpath 'Outputs'
addpath 'Structures'
addpath 'Tools'
addpath 'Structures/Materials'

% ----------------------------- Get time etc -----------------------------


cl = clock;                         % Get data and time for folder name

global R                            % Coefficient used in mapping
global S                            % Coefficient used in mapping
global s_par                        % Free space wavenumber
global p_par                        % Layer thickness array
global pola                         % Mode polarisation
global tol

tol = 1e-5;

% --------------------------- Section 1:  Read Input ---------------------

% Read initial parameters
in_par  =   input_parameters; 

% Read structure name
name    =   in_par{1};

% Read other parameters
par     =   in_par{2};

% Read initial photon parameters
lambda  =   par(1);                 % Free space wavelength (m)
omega   =   par(2);                 % Free space angular frequency (rad/s)
freq    =   par(3);                 % Free space frequency (1/s)
k0      =   par(4);                 % Free space wavevector (1/m)
E       =   par(5);                 % Free space photon energy (eV)

% Tests to run
m_in    =   par(6);             % Perform initial scan for modes (Yes/No)
rt_in   =   par(7);             % Perform initial reflection calculation (Yes/No)
m_scan  =   par(8);             % Track modes (number of runs)
rt_scan =   par(9);             % Track reflection trans (number of runs)

% Save parameters
prec    =   par(10);            % Number of significant figures to print
save_i  =   par(11);            % Save initial results to file (Yes/No)

% Modal analysis specific parameters

% Contour parameters
contour =   par(12:15);

% Other parameters
pola    =   par(16);                % Mode polarisation (TE or TM)
nrs     =   par(17);                % Number of Newton Raphson steps
vg      =   par(19);                % Calculate group velocity?
meth    =   par(20);                % Method

% Reflection transmission specific parameters
theta_i =   par(18);                % Initial angle of incidence (degrees)

% ------------------------ Section 2: Calculations ------------------------

p_par       =   [lambda omega freq k0 E];
s_par       =   struc(name,p_par);          % Structure parameter calculation

% Evaluate modes

if m_in == 1
    
    % Calculate structure parameters
    n_struc     =   s_par(3,:);
    nc          =   n_struc(1);
    ns          =   n_struc(end);
    R           =   ((ns^2) - (nc^2))/4;       % Coefficient used in mapping
    S           =   ((ns^2) + (nc^2))/2;        % Coefficient used in mapping
    
    % Find roots
    roots       =   find_roots(contour,nrs);
    
    clc
    fprintf('\nStage 1 Isolating roots  : Complete');
    fprintf('\nStage 2 Evaluating roots : Complete');
    
    % Save roots and initial conditions
    if save_i == 1;
        write_roots(cl,roots,name,prec,R,S,k0);
        fprintf('\n\nRoots saved\n');
    end
    
end

% Calculate reflection transmission coefficients if required

if rt_in == 1
    
    eps_struc   =   s_par(1,:);
    mu_struc    =   s_par(2,:);
    n_struc     =   s_par(3,:);
    z_struc     =   s_par(4,:);
    d_struc     =   s_par(5,3:end);
    
    co  =   coefficients(n_struc,z_struc,d_struc,theta_i,k0);
    
    R_par       =   co(1,2);
    T_par       =   co(1,1);
    A_par       =   1-(R_par+T_par);
    R_per       =   co(2,2);
    T_per       =   co(2,1);
    A_per       =   1-(R_per+T_per);
    
    ref         =   [R_par T_par A_par R_per T_per A_per];
    
    % Save roots and initial conditions
    if save_i == 1;
        write_r(cl,ref,name,prec,theta_i)
        fprintf('Reflection and transmission results saved');
    end
end



% Perform mode scans

if m_scan ~= 0
    if m_in ~= 0
        if save_i == 1
            ind = plot_roots(0,contour);
        else
            ind = plot_roots(roots,contour);
        end
    else
        ind = plot_roots(0,contour);
    end
    run = m_scan;    

    % Read roots from file
    FolderName  =   'tmp';
    FileName    =   'Modes.txt';
    FileName    =   sprintf('%s/%s',FolderName,FileName);
    
    FileID      =   fopen(FileName, 'r');
    loop = 0;
    
    while 1
        loop = loop+1;
        tline = fgetl(FileID);
        if ~ischar(tline),   break,   end
        
        if loop == 6
            num = sscanf(tline(16:end),'%g');
            a = zeros(num,7); 
        end
        
        if loop > 8
            a(loop-8,:) = sscanf(tline,'%g');
        end
    end
    roots = a;
    z_mat = complex(roots(ind,2),roots(ind,3));
    
      
    for loop = 1:run
        out_mat = m_tracker(in_par,run,loop,z_mat,cl,ind,vg,meth);
    end
    
    
end

% Perform reflection transmission scans

if rt_scan ~= 0    
    for loop = 1:rt_scan     
        R_par_mat = rt_tracker(in_par,rt_scan,loop,cl);
    end    
end


