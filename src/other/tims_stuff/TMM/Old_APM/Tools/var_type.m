function ind = var_type(var_name)


% Create variable strings

% Photon varaibles

str1 = {'lambda','wavelength'};
str2 = {'omega','angular'};
str3 = {'frequency','freq'};
str4 = {'k0','wavenumber','wave','wavevector'};
str5 = {'E','energy','photon'};

% Struc variables

str6 = {'permittivity','eps'};
str7 = {'permeabillity','mu'};
str8 = {'thickness','d'};

% Other

str9 = {'theta_i','incident','angle'};
str10 = {'beta'};


% Check variable name against above examples

var_name    =   [var_name ' '];
space       =   findstr(var_name,' ');
var_name    =   var_name(1:space(1));
var_name_n  =   deblank(var_name);

var(1) = max(strcmpi(var_name_n,str1));
var(2) = max(strcmpi(var_name_n,str2));
var(3) = max(strcmpi(var_name_n,str3));
var(4) = max(strcmpi(var_name_n,str4));
var(5) = max(strcmpi(var_name_n,str5));
var(6) = max(strcmpi(var_name_n,str6));
var(7) = max(strcmpi(var_name_n,str7));
var(8) = max(strcmpi(var_name_n,str8));
var(9) = max(strcmpi(var_name_n,str9));
var(10)= max(strcmpi(var_name_n,str10));


ind = find(var == 1);
