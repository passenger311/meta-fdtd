function out = least_squares(x,y)

n = length(x);
ssxx = sum(x.^2)-n*mean(x)^2;
ssyy = sum(y.^2)-n*mean(y)^2;
ssxy = sum(x.*y)-n*mean(x)*mean(y);

s = sqrt((ssyy-((ssxy^2)/ssxx))/(n-2));

m = ssxy/ssxx;
c = mean(y)-m*mean(x);

erc = s*sqrt((1/n)+(mean(x)^2)/ssxx);
erm = s/sqrt(ssxx);

%plot(x,y,'b',x,m*x+c,'r');

out = [m erm c erc];