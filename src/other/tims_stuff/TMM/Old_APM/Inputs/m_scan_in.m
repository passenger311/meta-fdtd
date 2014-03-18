function out = m_scan_in(run)


if run == 1
    c   =   299792458;
    lambdap =   (2*pi*c);
    % Mode scan 1 variables
    
    num_var     =   1;                  % Number of varaiables (1 or 2);
    st_var      =   1;                  % Change order of variables (1 or 2 first);
    
    var1_name   =   'omega';   % First variable (Any photon or layer parameter, or incident angle)
    var1_min    =   0.0001e15;                  % Min value of variable 1
    var1_max    =   2e15;                 % Max value of variable 1
    var1_inc    =   'Linear';           % Increment type of var 1 (Linear, Logarithmic or Custom)
    lin_inc1    =   [];                 % Linear increment (define either inc or res)
    lin_res1    =   1000;              % Linear resolution (number of points)
    log_res1    =   200;                % Log increment (number of points per decade)
    custom_inc1 =   read_custom;                 % Custom increment (define array of values)
    
    var2_name   =   'freq';    % Second variable (Any photon or layer parameter, or incident angle)
    var2_min    =   2e9;                  % Min value of variable 2
    var2_max    =   5e9;               % Max value of variable 2
    var2_inc    =   'Linear';           % Increment type of var 2 (Linear, Logarithmic or Custom)
    lin_inc2    =   [];                 % Linear increment (define either inc or res)
    lin_res2    =   2000;           % Linear resolution (number of points)
    log_res2    =   100;                % Log increment (number of points per decade)
    custom_inc2 =   read_custom;                 % Custom increment (define array of values)
    

    
    
elseif run == 2
    
    % Mode scan 2 variables
    
    num_var     =   1;                  % Number of varaiables (1 or 2);
    st_var      =   2;                  % Change order of variables (1 or 2 first);
    
    var1_name   =   'Permittivity of layer 1';   % First variable (Any photon or layer parameter, or incident angle)
    var1_min    =   4^2;                  % Min value of variable 1
    var1_max    =   5^2;                 % Max value of variable 1
    var1_inc    =   'Linear';           % Increment type of var 1 (Linear, Logarithmic or Custom)
    lin_inc1    =   [];                 % Linear increment (define either inc or res)
    lin_res1    =   [10];              % Linear resolution (number of points)
    log_res1    =   100;                % Log increment (number of points per decade)
    custom_inc1 =   [];                 % Custom increment (define array of values)
    
    var2_name   =   'Thickness of layer 2';    % Second variable (Any photon or layer parameter, or incident angle)
    var2_min    =   0.1e-9;                  % Min value of variable 2
    var2_max    =   1000e-9;               % Max value of variable 2
    var2_inc    =   'Log';           % Increment type of var 2 (Linear, Logarithmic or Custom)
    lin_inc2    =   [];                 % Linear increment (define either inc or res)
    lin_res2    =   lin_res1;           % Linear resolution (number of points)
    log_res2    =   100;                % Log increment (number of points per decade)
    custom_inc2 =   [];                 % Custom increment (define array of values)
    
elseif run == 3
    
    % Mode scan 3 variables
    
    num_var     =   2;                  % Number of varaiables (1 or 2);
    st_var      =   1;                  % Change order of variables (1 or 2 first);
    
    var1_name   =   'Incident angle';   % First variable (Any photon or layer parameter, or incident angle)
    var1_min    =   45;                 % Min value of variable 1
    var1_max    =   75;                 % Max value of variable 1
    var1_inc    =   'Linear';           % Increment type of var 1 (Linear, Logarithmic or Custom)
    lin_inc1    =   [];                 % Linear increment (define either inc or res)
    lin_res1    =   [300];              % Linear resolution (number of points)
    log_res1    =   100;                % Log increment (number of points per decade)
    custom_inc1 =   [];                 % Custom increment (define array of values)
    
    var2_name   =   'Photon energy';    % Second variable (Any photon or layer parameter, or incident angle)
    var2_min    =   3;                  % Min value of variable 2
    var2_max    =   4.25;               % Max value of variable 2
    var2_inc    =   'Linear';           % Increment type of var 2 (Linear, Logarithmic or Custom)
    lin_inc2    =   [];                 % Linear increment (define either inc or res)
    lin_res2    =   lin_res1;           % Linear resolution (number of points)
    log_res2    =   100;                % Log increment (number of points per decade)
    custom_inc2 =   [];                 % Custom increment (define array of values)
    
    % etc
    
end

var1    =   {var1_name var1_min var1_max var1_inc lin_inc1 lin_res1 log_res1 custom_inc1};
var2    =   {var2_name var2_min var2_max var2_inc lin_inc2 lin_res2 log_res2 custom_inc2};

if st_var == 2
    tmp     =   var2;
    var2    =   var1;
    var1    =   tmp;
end

out     =   {num_var st_var var1 var2};