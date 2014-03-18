function R_per_mat = rt_tracker(in_par,run_n,run,cl)
tic
% --------------------------- Section 1: Input ----------------------------
global pola

% Read structure name
name    =   in_par{1};

% Read other parameters
par     =   in_par{2};

% Read initial photon parameters
lambda  =   par(1);                % Free space wavelength (m)
omega   =   par(2);                 % Free space angular frequency (rad/s)
freq    =   par(3);                 % Free space frequency (1/s)
k0      =   par(4);                 % Free space wavevector (1/m)
E       =   par(5);                 % Free space photon energy (eV)

% Reflection transmission specific parameters
theta_i =   par(18);                % Initial angle of incidence (degrees)
prec    =   par(10);            % Number of significant figures to print

% Read variable parameters

var_in      =   rt_scan_in(run);

num_var     =   var_in{1};                  % Number of varaiables (1 or 2);


var1        =   var_in{3};

var1_name   =   var1{1};            % First variable (Any photon or layer parameter, or incident angle)
var1_min    =   var1{2};            % Min value of variable 1
var1_max    =   var1{3};            % Max value of variable 1
var1_inc    =   var1{4};            % Increment type of var 1 (Linear, Logarithmic or Custom)
lin_inc1    =   var1{5};            % Linear increment (define either inc or res)
lin_res1    =   var1{6};            % Linear resolution (number of points)
log_res1    =   var1{7};            % Log increment (number of points per decade)
custom_inc1 =   var1{8};            % Custom increment (define array of values)

var2        =   var_in{4};

var2_name   =   var2{1};            % Second variable (Any photon or layer parameter, or incident angle)
var2_min    =   var2{2};            % Min value of variable 2
var2_max    =   var2{3};            % Max value of variable 2
var2_inc    =   var2{4};            % Increment type of var 2 (Linear, Logarithmic or Custom)
lin_inc2    =   var2{5};            % Linear increment (define either inc or res)
lin_res2    =   var2{6};            % Linear resolution (number of points)
log_res2    =   var2{7};            % Log increment (number of points per decade)
custom_inc2 =   var2{8};            % Custom increment (define array of values)





% ------------------------ Section 2: Calculations ------------------------
c   =   299792458;

% Evaluate increment type

var1_inc    =   inc_type(var1_inc);
var2_inc    =   inc_type(var2_inc);


% Evaluate variable name type

ind1    =   var_type(var1_name);
ind2    =   var_type(var2_name);


% Calculate other photon parameters

p_par   =   [lambda,omega,freq,k0,E]; %    Photon parameters

% Calculate loop start and end values and increments

s1  =   scale(var1_min,var1_max,var1_inc,lin_inc1,lin_res1,log_res1,custom_inc1,2);
s2  =   scale(var2_min,var2_max,var2_inc,lin_inc2,lin_res2,log_res2,custom_inc2,num_var);

% Calculate structure parameters

s_par   =   struc(name,p_par);          % Structure parameter calculation


% Initialize loop1 counter

n1      =   0;


% Initialize results arrays

x_size      =   ceil((s1(2)-s1(1))/s1(3))+1;
y_size      =   ceil((s2(2)-s2(1))/s2(3))+1;

R_par       =   zeros(1,x_size);
T_par       =   zeros(1,x_size);
A_par       =   zeros(1,x_size);
R_per       =   zeros(1,x_size);
T_per       =   zeros(1,x_size);
A_per       =   zeros(1,x_size);

R_par_mat   =   zeros(y_size,x_size);
R_per_mat   =   zeros(y_size,x_size);
T_par_mat   =   zeros(y_size,x_size);
T_per_mat   =   zeros(y_size,x_size);
A_par_mat   =   zeros(y_size,x_size);
A_per_mat   =   zeros(y_size,x_size);

x_var       =   zeros(1,x_size);
y_var       =   zeros(1,y_size);


% Parameters used for caculating percentage of calculation complete

tot         =   x_size*y_size;
per_old     =   0;
sym         =   '%';
b = 0;

% Loop over variables and calculate coefficients

for loop1 = s2(1):s2(3):s2(2)
    
    if num_var == 2
        
        % Change order of variables if required and update parameters
        
        
        var1    =   var(loop1,var2_min,var2_inc,log_res2,custom_inc2,var2_name,ind2,s_par,p_par,theta_i,name,b);
        
        
        % Update other parameters
        
        theta_i =   var1(1,1);
        k0      =   var1(1,2);
        y       =   var1(2,1);
        b       =   var1(2,2);
        s_par   =   var1(3:end,:);
        p_par   =   photon_par([],[],[],k0,[]);
        
    else
        y       =   1;
    end
    
    n1      =   n1+1;
    n2      =   0;
    
    for loop2 = s1(1):s1(3):s1(2)
        
        n2      =   n2+1;
        
        % Change order of variables if required and update parameters
        
        
        var2    =   var(loop2,var1_min,var1_inc,log_res1,custom_inc1,var1_name,ind1,s_par,p_par,theta_i,name,b);
        
        
        % Update other parameters
        
        theta_i =   var2(1,1);
        k0      =   var2(1,2);
        omega   =   k0*c;
        x       =   var2(2,1);
        b       =   var2(2,2);
        s_par   =   var2(3:end,:);
        
        
        
        
        n_struc     =   s_par(3,:);
        z_struc     =   s_par(4,:);
        d_struc     =   s_par(5,3:end);
        
        
        eps_val   =   s_par(1,:);
        mu_val    =   s_par(2,:);
        
        d_val     =   s_par(5,3:end);
        
        nc          =   n_struc(1);
        ns          =   n_struc(end);
        R           =   ((ns^2) - (nc^2))/4;        % Coefficient used in mapping
        S           =   ((ns^2) + (nc^2))/2;        % Coefficient used in mapping
        
        if pola == 1
            m_val = s_par(2,:);
        elseif pola == 0;
            m_val = s_par(1,:);
        end
        
        
        
        % Calculate reflection and transmission coefficients
        
        co  =   coefficients(n_struc,z_struc,d_struc,theta_i,k0,b);
        %sol =   eigenfunction2(eps_val,d_val,m_val,mu_val,b,k0);
        
        R_par(n2)   =   co(1,2);
        T_par(n2)   =   co(1,1);
        A_par(n2)   =   1-(co(1,2)+co(1,1));
        R_per(n2)   =   co(2,2);
        T_per(n2)   =   co(2,1);
        A_per(n2)   =   1-(co(2,2)+co(2,1));
        
        x_var(n2)   =   x;
        
        
        % Calculate and print percentage of computation completed
        
        per         =   (round(100*(n1*n2)/tot))/1;
        
        if per > per_old+1
            clc
            fprintf('\nCalculation is %0.0f %s complete\n\n',per,sym)
            per_old = per_old+1;
        end
        
    end
    
    R_par_mat(n1,:) =   R_par;
    R_per_mat(n1,:) =   R_per;
    T_par_mat(n1,:) =   T_par;
    T_per_mat(n1,:) =   T_per;
    A_par_mat(n1,:) =   A_par;
    A_per_mat(n1,:) =   A_per;
    
    y_var(n1)   =   y;
    
end

% Trim arrays



R_par_mat   =   R_par_mat(1:n1,1:n2);
R_per_mat   =   R_per_mat(1:n1,1:n2);
T_par_mat   =   T_par_mat(1:n1,1:n2);
T_per_mat   =   T_per_mat(1:n1,1:n2);
A_par_mat   =   A_par_mat(1:n1,1:n2);
A_per_mat   =   A_per_mat(1:n1,1:n2);
n_struc     =   s_par(3,:);
n_struc(1)
x_var       =   x_var(1:n2);
y_var       =   y_var(1:n1);

if num_var == 2
    k0_mat      =   n_struc(1).*sin(x_var.*(pi/180));
    k0_mat      =   k0_mat'*y_var/(299792458*1e6);
    y_mat       =   ones(1,length(x_var));
    y_mat       =   y_mat'*y_var/1e15;
    for loop = 1:6
        filename    =   sprintf('Complex frequency scan outputs %g.txt',loop)
        fileID  =   fopen(filename,'r');
        loop2    =   0;
        data1   =   [];
        while 1
            loop2    =   loop2+1;
            tline   =   fgetl(fileID);
            
            if ~ischar(tline),   break, end
            
            a       =   textscan(tline,'%f');
            data1(loop2,:)   =   a{1}';
            
            
        end
        data_out{loop} = data1(:,1:2);
    end

    scl     =   'linear';
    figure(1)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(R_par_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Reflection (TM)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
        %plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])
%         figure(11)
%     set(gcf, 'Renderer', 'zbuffer');
%     surf(x_var,y_var,real(R_par_mat))
%     xlabel('k (pi/nd)','FontSize',16)
%     ylabel('k0 (pi/nd)','FontSize',16)
%     title('Reflection (TM)','FontSize',16)
%     shading interp
%     set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])

    
        figure(2)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(T_par_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Transmission (TM)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
        %plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])
%     figure(2)
%     set(gcf, 'Renderer', 'zbuffer');
%     surf(k0_mat,y_mat,real(T_par_mat'))
%     xlabel('k (pi/nd)','FontSize',16)
%     ylabel('k0 (pi/nd)','FontSize',16)
%     title('Transmission (TM)','FontSize',16)
%     shading interp
%     set(gca,'XScale',scl)
%     foo = get(gca,'dataaspectratio');
%     
%     set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])

    figure(3)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(A_par_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Absorption (TM)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
        %plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])

    figure(4)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(R_per_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Reflection (TE)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
        %plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])


    figure(5)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(T_per_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Transmission (TE)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
       % plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])

    figure(6)
    set(gcf, 'Renderer', 'zbuffer');
    pcolor(k0_mat,y_mat,real(A_per_mat'))
    xlabel('Longitudinal propagation constant, \beta (x10^6m^{-1})','FontSize',16)
    ylabel('Angular frequency, \omega (x10^{12} rad s^{-1})','FontSize',16)
    title('Absorption (TE)','FontSize',16)
    shading interp
    set(gca,'XScale',scl)
    %foo = get(gca,'dataaspectratio');
    %set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
    cl  =   {'k','g','r','b',[1 0.5 0.2],[0.5,0,0.5]};
    for loop = 1:6
        hold on
        dat1    =   data_out{loop};
        x   =   dat1(:,1)/1e6;
        y   =   dat1(:,2)/1e15;
       % plot(x,y,'Linewidth',2,'Color',cl{loop})
    end
    hold off
    %axis([0 100 0 800])

%     figure(6)
%     set(gcf, 'Renderer', 'zbuffer');
%     surf(k0_mat,y_mat,abs(real(A_par_mat')))
%     shading interp
%     set(gca,'XScale',scl)
%     xlabel('k (pi/nd)','FontSize',16)
%     ylabel('k0 (pi/nd)','FontSize',16)
%     title('Dispersion equation, F(w,b)','FontSize',16)
%     foo = get(gca,'dataaspectratio');
%     
%     set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)])
 %   view(2)
else
    figure(1)
    plot(x_var*1e-12,real(R_par_mat))
    figure(2)
    plot(x_var*1e-12,real(T_par_mat))
    figure(3)
    plot(x_var*1e-12,real(A_par_mat))
    figure(4)
    plot(x_var*1e-12,real(R_per_mat))
    figure(5)
    plot(x_var*1e-12,real(T_per_mat))
    figure(6)
    plot(x_var*1e-12,real(A_per_mat))
end


out_mat     =   {R_par_mat,R_per_mat,T_par_mat,T_per_mat,x_var,y_var};
output      =   {out_mat,prec,cl,name,run_n,run,k0,var1_name,var2_name,num_var};
%sav_rt(output);



t = toc;
hours = floor(t/3600);
mins  = floor(t/60)-(hours*60);
sec   = t - (mins*60)-(hours*3600);
clc
fprintf(1,'\nTime taken to run: %.0f hours %.0f mins %.1f secs\n\n',hours,mins,sec)

