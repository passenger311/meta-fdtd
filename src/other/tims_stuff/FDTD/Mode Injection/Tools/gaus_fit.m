function F = gauss_fit(p)

global x
global y

a = p(1);
w = p(2);
phi = p(3);
b = p(4);
c = p(5);

width = max(size(x));
offset = round(width/5);
r_width = width-(2*offset);
int = round(r_width/5);

for loop = 1:5
    ind = offset;
    x_mat(loop) = x(ind);
    y_mat(loop) = y(ind);
    offset = offset+int;
end


f1 = (a*sin((w*x_mat(1))+phi)*exp(-((x_mat(1)-b)^2)/(2*(c^2))))-y_mat(1);
f2 = (a*sin((w*x_mat(2))+phi)*exp(-((x_mat(2)-b)^2)/(2*(c^2))))-y_mat(2);
f3 = (a*sin((w*x_mat(3))+phi)*exp(-((x_mat(3)-b)^2)/(2*(c^2))))-y_mat(3);
f4 = (a*sin((w*x_mat(4))+phi)*exp(-((x_mat(4)-b)^2)/(2*(c^2))))-y_mat(4);
f5 = (a*sin((w*x_mat(5))+phi)*exp(-((x_mat(5)-b)^2)/(2*(c^2))))-y_mat(5);

F = [f1;f2;f3;f4;f5];