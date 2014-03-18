function out = pow3(num)

    pow     =   floor(log(num)/log(1000));
    num     =   num*10^(pow*-3);
    pref    =   'yzafpnum kMGTPEZY';
    out     =   {num,pref(pow+9),pow*-3};
end
    
   
    
    
    