close all
clear all
pause(1)

c           =   299792458;
res         =   10;
lambda      =   1.55e-6;
eps_max     =   11.68;
dx          =   lambda/res/sqrt(eps_max);
cour        =   0.7;
dt          =   cour*dx/c;

Filename1   =   'ebal_EnI2.0.gpl';
Filename2   =   'ebal_divs1.0.gpl';
Filename3   =   'ebal_divs3.0.gpl';


fileID  =   fopen(Filename1,'r');

for loop = 1:13
    tline   =   fgetl(fileID);
    if loop == 11
        a   =   textscan(tline,'%s');
        params(1)   =   str2double(a{1}(2));
        params(2)   =   str2double(a{1}(3));
    end
end

loop    =   0;
while 1
    loop    =   loop+1;
    tline   =   fgetl(fileID);
    
    if ~ischar(tline),   break, end
    
    a       =   textscan(tline,'%f');
    data1(loop,:)   =   a{1}';
    
end
fclose(fileID);
% figure(1)
% plot(data1(:,1))
% hold off
% figure(2)
x   =   dt.*(1:length(data1(:,1)));
y   =   log(data1(:,1));
% plot(log(data1(:,1)))
m_y = mean((y(2:end)-y(1:end-1))./dt);
c_y = mean(y'-(m_y.*x));
s_params = [m_y, c_y];
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
estimate = fminsearch(@linfit,s_params,options,x,y');
m_c   =   estimate(1);
c_c   =   estimate(2);


fileID  =   fopen(Filename2,'r');
for loop = 1:13
    tline   =   fgetl(fileID);
end

loop    =   0;
while 1
    loop    =   loop+1;
    tline   =   fgetl(fileID);
    
    if ~ischar(tline),   break, end
    
    a       =   textscan(tline,'%f');
    data2(loop,:)   =   a{1}';
    
end

fclose(fileID);

% figure(3)
% plot(data2(:,1))
% figure(4)
% plot(data2(:,2))
% figure(5)
% plot(data2(:,3))
% figure(6)
% plot(data2(:,3)./data1(:,1))
wat     =   data2(:,3)*((c^3)/dx);
joul    =   data1(:,1)*(c^2);
lft     =   wat./joul;
mlft    =   mean(lft);


fileID  =   fopen(Filename3,'r');
for loop = 1:13
    tline   =   fgetl(fileID);
end

loop    =   0;
while 1
    loop    =   loop+1;
    tline   =   fgetl(fileID);
    
    if ~ischar(tline),   break, end
    
    a       =   textscan(tline,'%f');
    data3(loop,:)   =   a{1}';
    
end

fclose(fileID);

% figure(7)
% plot(data3(:,1))
% figure(8)
% plot(data3(:,2))
% figure(9)
% plot(data3(:,3))
% figure(10)
% plot(data3(:,3)./data1(:,1))
wat     =   data3(:,3)*((c^3)/dx);
joul    =   data1(:,1)*(c^2);
rgt     =   wat./joul;
mrgt    =   mean(rgt);
diff    =   mrgt-mlft;
vg      =   (mrgt-mlft)*dx/c;
vg      =   vg*(params(2)-params(1));


filename    =   'Energy_vel.txt';
fileID      =   fopen(filename,'w');

fprintf(fileID,'Box width               =   %g grid cells\n',params(2)-params(1));
fprintf(fileID,'Energy flow to right    =   %g s-1\n',mrgt);
fprintf(fileID,'Energy flow to left     =   %g s-1\n',mlft);
fprintf(fileID,'Difference              =   %g s-1\n',diff);
fprintf(fileID,'Vg(maybe)               =   %g c\n',vg);

fclose(fileID);
%quit;
