function par = var(step,var_min,var_inc,log_inc,custom_inc,var_name,ind,s_par,p_par,theta_i,name,b)

c       =   299792458;
if var_inc == 1
    par = step;
elseif var_inc == 2
    par = (10^((step-1)/log_inc))*var_min;
elseif var_inc == 3
    par = custom_inc(step);
end


if size(ind) > 1 | isempty(ind) == 1 
    fprintf(['Error: variable name not recognised, please see function '...
        ,'"var.m" for list of variable names\n\n']);
    return
elseif ind == 1
    lambda  =   par;
    p_par   =   photon_par(lambda,[],[],[],[]);
    s_par   =   struc(name,p_par);
   
elseif ind == 2
    omega   =   par;
    p_par   =   photon_par([],omega,[],[],[]);
    s_par   =   struc(name,p_par);
elseif ind == 3
    freq    =   par;
    p_par   =   photon_par([],[],freq,[],[]);
    s_par   =   struc(name,p_par);
elseif ind == 4
    k0      =   par;
    p_par   =   photon_par([],[],[],k0,[]);
    s_par   =   struc(name,p_par);
elseif ind == 5
    E       =   par;
    p_par   =   photon_par([],[],[],[],E);
    s_par   =   struc(name,p_par);
elseif ind == 6
    nums    =   regexp(var_name,'\d+','match');
    layer   =   str2double(nums);
          
    if isempty(findstr(var_name,'cover')) == 0
        layer   =   [layer 1];
    end

    if isempty(findstr(var_name,'substrate')) == 0
        layer   =   [layer size(s_par,2)];
    end
         
    s_par(1,layer)  =   par;
    s_par(3,:)      =   sqrt(s_par(1,:).*s_par(2,:));
    s_par(4,:)      =   sqrt(s_par(2,:)./s_par(1,:));
elseif ind == 7
    nums    =   regexp(var_name,'\d+','match');
    layer   =   str2double(nums);
    
    if isempty(findstr(var_name,'cover')) == 0
        layer   =   [layer 1];
    end

    if isempty(findstr(var_name,'substrate')) == 0
        layer   =   [layer size(s_par,2)];
    end
    s_par(2,layer)  =   par;
    s_par(3,:)      =   sqrt(s_par(1,:).*s_par(2,:));
    s_par(4,:)      =   sqrt(s_par(2,:)./s_par(1,:)); 
elseif ind == 8
    nums    =   regexp(var_name,'\d+','match');
    layer   =   str2double(nums);
    s_par(5,layer+1)    =   par;
elseif ind == 9
    theta_i = par;
elseif ind == 10

    b   =   par;
end

par1        =   zeros(2,size(s_par,2));
par1(1,1)   =   theta_i;
par1(1,2)   =   p_par(4);
par1(2,1)   =   par;
par1(2,2)   =   b;

par     =   [par1;s_par];
    
