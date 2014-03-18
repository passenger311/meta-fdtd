close all
clear all
pause(1)
clc

c           =   299792458;
res         =   10;
lambda      =   1.55e-6;
eps_max     =   11.68;
dx          =   0.01e-6;
cour        =   0.7;
dt          =   cour*dx/c;

Filename{1} =   'ebal_energy.0.gpl';
Filename{2} =   'ebal_energy_int.0.gpl';
Filename{3} =   'ebal_flux.0.gpl';

for loop    =   1:3
    
    fileID  =   fopen(Filename{loop},'r');
    loop2   =   0;
    data    =   [];
    for l   =   1:13
        tline   =   fgetl(fileID);
        if l == 8
            inf     =   textscan(tline,'%s');
            mult    =   str2double(inf{1}{4});
        end
    end
    while 1
        loop2   =   loop2+1;
        tline   =   fgetl(fileID);
        
        if ~ischar(tline),   break, end
        
        a               =   textscan(tline,'%f');
        data(loop2,:)   =   a{1}';
        
    end
    data(1:3,:)     =   [];
    data_out{loop}  =   data;
    fclose(fileID);
end

En  =   data_out{1};
EnI =   data_out{2};
Ds  =   data_out{3};
tot_enI     =   EnI(:,1);
flux    =   Ds(:,2);
time    =   (1:length(En))*dt*mult;
cl  =   {'b','r','k','g',[1 0.5 0.2]};
figure(1)
for loop = 1:5
    plot(time*1e15,(En(:,loop))*(c^3)/dx,'color',cl{loop});
    hold on
end
hold off
title('En');
xlabel('Time (fs)')
ylabel('Energy flux (Js^{-1})')
set(gcf,'color','white');
legend('dudt','ds','jebas','khabs','res');
filename    =   'Results/En.eps';
print(gcf,'-depsc','-r200','-painters',filename);
figure(2)
for loop = 1:5
    plot(time*1e15,(EnI(:,loop))*c^2,'color',cl{loop});
    hold on
end
hold off
title('EnI');
xlabel('Time (fs)')
ylabel('Energy (J)')
set(gcf,'color','white');
legend('sumdudt','sumds','sumje','sumkh','sumres');
filename    =   'Results/EnI.eps';
print(gcf,'-depsc','-r200','-painters',filename);
figure(3)
for loop = 1:4
    plot(time*1e15,(Ds(:,loop))*(c^3)/dx,'color',cl{loop});
    hold on
end
hold off
title('Ds');
xlabel('Time (fs)')
ylabel('Energy flux (Js^{-1})')
set(gcf,'color','white');
legend('ds','dsx','dsy','dsz');
filename    =   'Results/Ds.eps';
print(gcf,'-depsc','-r200','-painters',filename);
figure(4)
plot(time*1e15,(Ds(:,2)./EnI(:,1))*c/dx);
set(gcf,'color','white');
title('Ds(x)/EnI, Radiative (s^{-1})')
xlabel('Time (fs)')
ylabel('Energy Radiated (s^{-1})')
filename    =   'Results/Radiative_x.eps';
print(gcf,'-depsc','-r200','-painters',filename);
figure(5)
plot(time*1e15,(Ds(:,3)./EnI(:,1))*c/dx);
set(gcf,'color','white');
title('Ds(y)/EnI, Radiative (s^{-1})')
xlabel('Time (fs)')
ylabel('Energy Radiated (s^{-1})')
filename    =   'Results/Radiative_y.eps';
print(gcf,'-depsc','-r200','-painters',filename);
figure(6)
plot(time*1e15,(En(:,3)./EnI(:,1))*c/dx);
set(gcf,'color','white');
title('En/EnI, Dispersion (s^{-1})')
xlabel('Time (fs)')
ylabel('Energy dispersion (s^{-1})')
filename    =   'Results/Dispersion.eps';
print(gcf,'-depsc','-r200','-painters',filename);

leakyx   =   (Ds(:,2)./EnI(:,1));
leakyy  =   (Ds(:,3)./EnI(:,1));
disp    =   (En(:,3)./EnI(:,1));

filename    =   'Results/Loss.txt';
fileID  =   fopen(filename,'w');
gamma_leakyx     =   mean(leakyx)*c/dx
st_leakyx   =   std(leakyx)*c/dx
gamma_leakyy     =   mean(leakyy)*c/dx
st_leakyy   =   std(leakyy)*c/dx
gamma_disp      =   mean(disp)*c/dx
st_disp   =   std(disp)*c/dx
ratio   =   gamma_leakyx*100/gamma_disp

fprintf(fileID,'Radiative energy loss (x) = %g +/- %g s^-1\n',gamma_leakyx,st_leakyx);
fprintf(fileID,'Radiative energy loss (y) = %g +/- %g s^-1\n',gamma_leakyy,st_leakyy);
fprintf(fileID,'Dispersive energy loss = %g +/- %g s-1\n',gamma_disp,st_disp);
fprintf(fileID,'Loss ratio R/D = %g',ratio);
fclose(fileID);
quit;