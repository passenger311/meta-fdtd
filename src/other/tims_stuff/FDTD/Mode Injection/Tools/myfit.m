function sse = myfit(params,x,y)

a = params(1);
w = params(2);
phi = params(3);
b = params(4);
c = params(5);
%  a2 = params(6);
%  w2 = params(7);
%  phi2 = params(8);
%  b2 = params(9);
%  c2 = params(10);

fitted_curve = a.*sin((w.*x)+phi).*exp(-((x-b).^2)./(2*(c^2)));
%  fitted_curve2 = a2.*sin((w2.*x)+phi2).*exp(-((x-b2).^2)./(2*(c2^2)));
%  fitted_curve = fitted_curve1+fitted_curve2;

error_vector = fitted_curve-y;
sse = sum(error_vector.^2);