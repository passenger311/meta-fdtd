function sse = myfit2(params,x,y)

a = params(1);
b = params(2);
c = params(3);

fitted_curve = a.*exp(-((x-b).^2)./(2*(c^2)));

error_vector = fitted_curve-y;
sse = sum(error_vector.^2);