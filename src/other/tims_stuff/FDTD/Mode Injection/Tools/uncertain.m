function out = uncertain(num,err)
    pow     =   abs(floor(log(err)/log(10)));
    err     =   rnd(err,pow);
    num     =   rnd(num,pow);
    out     =   [num,err];
end