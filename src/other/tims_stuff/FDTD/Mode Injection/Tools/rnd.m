function out = rnd(num,idp)
    mult    =   10^idp;
    out     =   floor(num*mult+0.5)/mult;
end