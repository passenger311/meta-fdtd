function out = lin_fit(x,y)

% Calculate guess values for the gradient and intercept
guess.m   =   mean(y(2:end)-y(1:end-1));
guess.c   =   mean(y'-guess.m*x);

% Find best fit values using fminsearch function
s_params = [guess.m, guess.c];
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
estimate = fminsearch(@linfit,s_params,options,x,y');

best.m  =   estimate(1);
best.c  =   estimate(2);


% Calculate error using parallelogram method
u_fit = (best.m.*x)+best.c;
u_diff = y-u_fit';
a_new   =   max(abs(u_diff))/2;
ex = 0;
loop = 0;
while ex == 0
    loop = loop +1;
    ind =   find(u_diff>=-a_new & u_diff<=a_new);
    nm  =   length(ind)/length(y);
    if (nm >= 0.68 & nm <=0.69);
        ex = 1;
    elseif nm < 0.68
        a_new = a_new*1.001;
    elseif nm > 0.69
        a_new = a_new*0.999;
    end
    
    if loop == 5000
        break
    end
    
end

par_u   =   u_fit+a_new;
par_l   =   u_fit-a_new;
dx  =   x(end)-x(1);
dy  =   par_u(end)-par_l(1);
m_u =   dy/dx;
c_u =   par_l(1)-x(1)*m_u;
dy  =   par_l(end)-par_u(1);
m_l =   dy/dx;
c_l =   par_u(1)-x(1)*m_l;

error.m     =   abs(m_u-m_l)/2;
error.c     =   abs(c_u-c_l)/2;


out     =   {best.m,error.m,best.c,error.c};
    
end