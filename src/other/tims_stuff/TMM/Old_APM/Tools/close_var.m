function out = close_var(var,var_min,var_inc,log_inc,custom_inc,s)

num     =   ceil((abs(s(1)-s(2)))/s(3));
par_mat =   zeros(1,num);
loop    =   0;

for step = s(1):s(3):s(2)
    loop = loop+1;
    if var_inc == 1
        par = step;
    elseif var_inc == 2
        par = (10^((step-1)/log_inc))*var_min;
    elseif var_inc == 3
        par = custom_inc(step);
    end
    par_mat(loop) = par;
end

diff    =   abs(par_mat-var(1));

ind     =   find(diff==min(diff));
ind     =   ind(1);

out = [ind,par_mat(ind)];


