% ------------------------------------------------------------------------
%           Transfer Matrix Method (TMM) Mode Solver Demonstration
% 
% This program demonstrates some of the features of the TMM mode solver
% and associated waveguide diagnostic functions contained in the class
% definition 'Structure.m'. In the first section definition of material 
% layers is introduced before solving the dispersion equation for complex
% frequency and complex wavevector solutions. Following this initial 
% section reflection transmission and absorption calculations are shown.
% Finally a method for calculating the energy confinement factor, energy
% velocity and also dissipative and radiative contributions to the total
% loss rate are shown.
% ------------------------------------------------------------------------


%% 1. Initialising the structure 

% For reference the double %% marks the start of a new code section which
% can be evaluated seperately from the rest of the program by pressing
% crtl + enter

close all
clear all
clc

% If using the four-level gain model the following parameters have to be
% defined:

dens = 2e18;        % Carrier density (cm^-3)
etau = 40e-15;      % Emission decay rate
atau = 40e-15;      % Absorption decay rate
elam = 1.55e-6;     % Emission wavelength (m)
alam = 1.345e-6;    % Absorption wavelength (m)
epsh = 11.68;       % Host permittivity
ed = 0e-9;          % Emission 
ad = 0e-9;          % Absorption
inv = 0;            % Inverison of lasing transition


% Initialise structure

waveguide = Structure({...              % Define Layers
    'Metal(4,3.13e15,1.07e13)',...      % Dispersive material models such as Drude(eps host, wp, gamma)
    'Dielectric(11.68,1)',...         	% General Dielectric(eps,mu)
    sprintf('four_lvl(%g,%g,%g,%g,%g,%g,%g,%g,%g)',epsh,alam,elam,ad,ed,1/atau,1/etau,dens*1e6,inv),... % Four level gain
    'Si',...
    'Metal(4,3.13e15,1.07e13)',...
    });

waveguide.d = [10 270 10].*1e-9;     % Layer thicknesses (m)
waveguide.pola = 'TM';                       % Field polarisation (TE/TM)       
waveguide.bnd = 0;                           % Mode type (0=bound, 1=leaky(sup), 2=leaky(sub), 3=leaky(both))

%% 2. Plot Dispersion equation

close all
clc

% A useful first step when analysing a new waveguide is to plot the
% dispersion equation in order to visulaise the modal dispersion and see
% roughly how many modes there are and where abouts they lie on the
% dispersion plot.

% This is done using the function realdispersion

b_min = 0;      % Minimum value for beta range
b_max = 4e7;    % Maximum value for beta range
b_res = 200;    % Beta range resolution
w_min = 0;      % Same for angular frequency
w_max = 5e15;
w_res = 200;

disp = waveguide.realdispersion(b_min,b_max,b_res,w_min,w_max,w_res);
% The output is a cell array with the beta range, omega range, and the 2D
% dispersion function.

% Plot the dispersion function (I've found the log(abs(dispersion)) gives
% the best picture of the root locations (shown by dark lines)

h = figure(1);
set(h,'color','w')
pcolor(disp{1}/1e6,disp{2}/1e15,log(abs(disp{3})))
shading interp
colormap(hot(1024))
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');


%% 3. Calculate the complex-frequency dispersion

% From the rough plot of the modal dispersion calculated in section 2 a
% rough estimate of the modal positions was seen. Here the dispersion
% equation is solved to find the complex-frequency roots using a
% Newton-Raphson method requiring initial guess values. Once these roots
% are found the modal fields can be plotted to characterise the roots, the
% dispersion curve of selected roots can then be traced over a defined beta
% range, again using the Newton-Raphson method.
clc
% Step 1: find the roots

% Define range of initial guess values
b_in  = 10e6;       % Initial beta value
w_min = 0;          % Minimum value for omega range
w_max = 1.5e15;       % Maximum value for omega range
w_res = 200;        % Omega range resolution

% Find roots (output: column 1 root number, column 2 complex root value
tmp = waveguide.find_roots_w(w_min,w_max,w_res,b_in);
roots = tmp(:,2);

% Print roots to screen
fprintf('\n%g roots found:\n',length(roots));
for n = 1:length(roots)
    fprintf('\nRoot %g has freq = %g THz and loss = %g ps^-1',n,real(roots(n))/2/pi/1e12,...
        imag(roots(n))*2/1e12);
end
fprintf('\n');

% Ask for roots to plot
a = input('\n\nPlease select root to plot: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end

% Plot field profile
close all
for loop = 1:length(a)
    % Field_plotter input (angular frequency, beta, length into
    % superstrate, length into substrate, resolution, plot(0=off,1=on)
    field = waveguide.field_plotter(roots(a(loop)),b_in,200e-9,500e-9,500,1);
    % Output: {1} = 4 columns x_mat|z_field|y_field|x_field
    %         {2} = permittivity,{3} = permittivity dispersion, {4} =
    %         permeability {5} = permeability disperison
end

% Step 2: calculate dispersion curves:

% Ask for roots to track
a = input('\n\nPlease select root(s) to track: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end
close all

% Track roots
% Over a beta range containing the initial beta value where the root was
% found
b_min = 0;
b_max = 40e6;
b_res = 500;

for i = 1:length(a)
    fprintf(1,sprintf('\nTracking mode %g of %g\n',i,length(a)));
    mode = waveguide.dispersionW(roots(a(i)),b_in,b_min,b_max,b_res);
    b_mat(:,i) = mode(:,1);
    w_mat(:,i) = mode(:,3);
    imw_mat(:,i) = mode(:,4);
    loss_mat(:,i) = -2*mode(:,4);
    vg_mat(:,i) = real(mode(:,5));
end

% Plot results
cl = lines(length(a));
for i = 1:length(a)
    h1 = figure(1);
    set(h1,'color','w')
    
    subplot(2,2,1)
    hold on
    plot(b_mat(:,i)/1e6,w_mat(:,i)/1e15,'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})')
    hold off
    
    subplot(2,2,2)
    hold on
    plot(b_mat(:,i)/1e6,loss_mat(:,i)/1e12,'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Modal loss rate \gamma (x10^{12}s^{-1})')
    hold off
    
    subplot(2,2,3)
    if i > 1
        hold on     % issue with semilog plot and hold on
    end
    semilogy(b_mat(:,i)/1e6,abs(vg_mat(:,i)),'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Group velocity v_{g} (c)')
    hold off
end

% Save results to file
for i = 1:length(a)
    filename = sprintf('Complex_w_%s_Mode_%g.txt',waveguide.pola,a(i));
    M = [b_mat(:,i) w_mat(:,i) imw_mat(:,i) loss_mat(:,i) vg_mat(:,i)];
    dlmwrite(filename,M,'\t');
end

%% 4. Calculate the complex-beta dispersion

% From the rough plot of the modal dispersion calculated in section 2 a
% rough estimate of the modal positions was seen. Here the dispersion
% equation is solved to find the complex-beta roots using a
% Newton-Raphson method requiring initial guess values. Once these roots
% are found the modal fields can be plotted to characterise the roots, the
% dispersion curve of selected roots can then be traced over a defined
% omega range, again using the Newton-Raphson method.
clc
% Step 1: find the roots

% Define range of initial guess values
w_in  = 1.4e15;     % Initial angular frequency value
b_min = 0;          % Minimum value for beta range
b_max = 20e6;     % Maximum value for beta range
b_res = 200;        % Beta range resolution

% Find roots (output: column 1 root number, column 2 complex root value
tmp = waveguide.find_roots_k(b_min,b_max,b_res,w_in);
roots = tmp(:,2);

% Print roots to screen
fprintf('\n%g roots found:\n',length(roots));
for n = 1:length(roots)
    fprintf('\nRoot %g has beta = %g x10^6m^-1 and loss = %g cm^-1',n,real(roots(n))/1e6,...
        imag(roots(n))*2/1e2);
end
fprintf('\n');

% Ask for roots to plot
a = input('\n\nPlease select root(s) to plot: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end

% Plot field profile
close all
for loop = 1:length(a)
    % Field_plotter input (angular frequency, beta, length into
    % superstrate, length into substrate, resolution, plot(0=off,1=on)
    field = waveguide.field_plotter(w_in,roots(a(loop)),200e-9,500e-9,500,1);
    % Output: {1} = 4 columns x_mat|z_field|y_field|x_field
    %         {2} = permittivity,{3} = permittivity dispersion, {4} =
    %         permeability {5} = permeability disperison
end

% Step 2: calculate dispersion curves:

% Ask for roots to track
a = input('\n\nPlease select root(s) to track: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end
close all

% Track roots
% Over a omega range containing the initial omega value where the root was
% found
w_min = 0;
w_max = 2e15;
w_res = 500;

for i = 1:length(a)
    fprintf(1,sprintf('\nTracking mode %g of %g\n',i,length(a)));
    mode = waveguide.dispersionK(roots(a(i)),w_in,w_min,w_max,w_res);
    w_mat(:,i) = mode(:,1);
    b_mat(:,i) = mode(:,3);
    imb_mat(:,i) = mode(:,4);
    loss_mat(:,i) = 2*mode(:,4);
    vg_mat(:,i) = real(mode(:,5));
end

% Plot results
cl = lines(length(a));
for i = 1:length(a)
    h1 = figure(1);
    set(h1,'color','w')
    
    subplot(2,2,1)
    hold on
    plot(b_mat(:,i)/1e6,w_mat(:,i)/1e15,'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})')
    hold off
    
    subplot(2,2,2)
    hold on
    plot(loss_mat(:,i)/1e2,w_mat(:,i)/1e15,'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Modal loss rate \gamma (x10^{12}s^{-1})')
    hold off
    
    subplot(2,2,3)
    if i > 1
        hold on     % issue with semilog plot and hold on
    end
    plot((vg_mat(:,i)),w_mat(:,i)/1e15,'color',cl(i,:),'linewidth',2);
    xlabel('Propagation constant \beta (x10^{6}m^{-1})');
    ylabel('Group velocity v_{g} (c)')
    xlim([-1 1])
    hold off
end

% Save results to file
for i = 1:length(a)
    filename = sprintf('Complex_k_%s_Mode_%g.txt',waveguide.pola,a(i));
    M = [w_mat(:,i) b_mat(:,i) imb_mat(:,i) loss_mat(:,i) vg_mat(:,i)];
    dlmwrite(filename,M,'\t');
end

%% 5. Calculate reflection transmission absorption data

% 2D plot
w_min = 0;
w_max = 5e15;
w_res = 200;
b_min = 0;
b_max = 40e6;
b_res = 200;
data = waveguide.reflec2D(w_min,w_max,w_res,b_min,b_max,b_res,'Beta'); % Can also use angle
b_mat = data{1};
w_mat = data{2};
TM_R = real(data{3});
TM_T = real(data{4});
TM_A = real(data{5});
TE_R = real(data{6});
TE_T = real(data{7});
TE_A = real(data{8});


close all
h = figure(1);
set(h,'color','w')
subplot(2,2,1)
pcolor(b_mat/1e6,w_mat/1e15,TM_R)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TM Reflection')
subplot(2,2,2)
pcolor(b_mat/1e6,w_mat/1e15,TM_T)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TM Transmission')
subplot(2,2,3)
pcolor(b_mat/1e6,w_mat/1e15,TM_A)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TM Absorption')

h = figure(2);
set(h,'color','w')
subplot(2,2,1)
pcolor(b_mat/1e6,w_mat/1e15,TE_R)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TE Reflection')
subplot(2,2,2)
pcolor(b_mat/1e6,w_mat/1e15,TE_T)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TE Transmission')
subplot(2,2,3)
pcolor(b_mat/1e6,w_mat/1e15,TE_A)
shading interp
colormap(hot(1024))
caxis([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('Angular frequency \omega (x10^{15}rad\cdots^{-1})');
title('TE Absorption')

% 1D plot
waveguide.beta = 10e6;
w_min = 0;
w_max = 5e15;
w_res = 2000;
data = waveguide.reflec1D(w_min,w_max,w_res,'Omega');
w_mat = data(:,1);
TM_R = real(data(:,2));
TM_T = real(data(:,3));
TM_A = real(data(:,4));
TE_R = real(data(:,5));
TE_T = real(data(:,6));
TE_A = real(data(:,7));

h = figure(3);
set(h,'color','w')
subplot(2,1,1)
plot(w_mat/1e15,TM_R,'r','linewidth',2)
hold on
plot(w_mat/1e15,TM_T,'linewidth',2)
plot(w_mat/1e15,TM_A,'g','linewidth',2)
ylim([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('RTA');
title('TM')
subplot(2,1,2)
plot(w_mat/1e15,TE_R,'r','linewidth',2)
hold on
plot(w_mat/1e15,TE_T,'linewidth',2)
plot(w_mat/1e15,TE_A,'g','linewidth',2)
ylim([0 1])
xlabel('Propagation constant \beta (x10^{6}m^{-1})')
ylabel('RTA');
title('TE')

