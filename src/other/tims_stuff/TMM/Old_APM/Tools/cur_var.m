function par = cur_var(var_name,ind,s_par,phot)

lambda  =   phot(1);
omega   =   phot(2);
freq    =   phot(3);
k0      =   phot(4);
E       =   phot(5);

if size(ind) > 1 | isempty(ind) == 1 
    fprintf(['Error: variable name not recognised, please see function '...
        ,'"var.m" for list of variable names\n\n']);
    return
elseif ind == 1
    par     =   lambda;
elseif ind == 2
    par     =   omega;
elseif ind == 3
    par     =   freq;
elseif ind == 4
    par     =   k0;
elseif ind == 5
    par     =   E;
elseif ind == 6
    nums    =   regexp(var_name,'\d+','match');
    layer   =   (str2double(nums));        
    if isempty(findstr(var_name,'cover')) == 0
        layer   =   [layer 1];
    end
    if isempty(findstr(var_name,'substrate')) == 0
        layer   =   [layer size(s_par,2)];
    end         
    par             =   s_par(1,layer);
elseif ind == 7
    nums    =   regexp(var_name,'\d+','match');
    layer   =   (str2double(nums));    
    if isempty(findstr(var_name,'cover')) == 0
        layer   =   [layer 1];
    end

    if isempty(findstr(var_name,'substrate')) == 0
        layer   =   [layer size(s_par,2)];
    end
    par     =   s_par(2,layer);
elseif ind == 8
    nums    =   regexp(var_name,'\d+','match');
    layer   =   str2double(nums);
    par     =   s_par(5,layer+1);
elseif ind == 9
    fprintf('Error: Please do not select incident angle as a variable');
    return
end



    
