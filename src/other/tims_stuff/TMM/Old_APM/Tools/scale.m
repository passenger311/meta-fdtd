function s = scale(var_min,var_max,var_inc,lin_inc,lin_res,log_inc,custom_inc,num_var)


if num_var == 2
    if var_inc == 1
        st  =   var_min;
        en  =   var_max;
        if isempty(lin_inc) == 1
            inc     =   (en-st)/lin_res;
        elseif isempty(lin_res) == 1
            inc =   lin_inc;
        else
            fprintf('\nError: Please define only lin_inc or lin_res\n\n');
        end
    elseif var_inc == 2
        st  =   1;
        en  =   var_max/var_min;
        en  =   log10(en);
        en  =   (en*log_inc)+1;
        en  =   round(en);
        inc =   1;
    elseif var_inc == 3
        st  =   1;
        en  =   length(custom_inc);
        inc =   1;
    end
else
    st  = 1;
    en  = 1;
    inc = 1;
end

s   =   [st en inc];
