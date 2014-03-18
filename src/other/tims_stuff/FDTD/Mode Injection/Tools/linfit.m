function sse = linfit(params,x,y)

m = params(1);
c = params(2);


fitted_curve = (m.*x)+c;

error_vector = fitted_curve-y;
sse = sum(error_vector.^2);