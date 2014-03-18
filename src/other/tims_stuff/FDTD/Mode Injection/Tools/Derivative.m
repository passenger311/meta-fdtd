close all
clear all
clc

c   =   299792458;

filename    =   'Evel_fit.txt';
fileID  =   fopen(filename,'r');

[A,COUNTA] = textscan(fileID,'%n',1e5,'Delimiter','\n','CommentStyle','#');

A   =   A{1};

t   =   A(1:3:end);
y   =   A(2:3:end)*1e-8;


d_y =   diff(y)./diff(t);
vg  =   d_y/c;
vg_m    =   mean(vg);
vg_sd   =   std(vg);

a   =   sprintf('The average group velocity is Vg = %g +/- %g\n',vg_m,vg_sd);
fprintf(a)
vg_r    =   uncertain(vg_m,vg_sd);
a   =   sprintf('The average group velocity is Vg = %g +/- %g\n',vg_r(1),vg_r(2));
fprintf(a)


time    =   pow3(max(t));
yn      =   pow3(max(y));
vgn     =   pow3(max(vg));

figure(3)
subplot(3,1,1)
plot(t*(10^time{3}),y*(10^yn{3}))
xlabel(sprintf('Time (%ss)',time{2}));
ylabel(sprintf('Centre position, (%sm)',yn{2}))
v   =   axis;
axis([min(t*(10^time{3})),max(t*(10^time{3})),v(3),v(4)])
subplot(3,1,2)
plot(t(1:end-1)*(10^time{3}),vg*(10^vgn{3}))
xlabel(sprintf('Time (%ss)',time{2}));
ylabel(sprintf('Energy velocity, x10^%s (C)',num2str(-vgn{3})))
v   =   axis;
axis([min(t(1:end-1)*(10^time{3})),max(t(1:end-1)*(10^time{3})),v(3),v(4)])
subplot(3,1,3)
plot(t(2:end)*(10^time{3}),vg*(10^vgn{3}))
xlabel(sprintf('Time (%ss)',time{2}));
ylabel(sprintf('Energy velocity, x10^%s (C)',num2str(-vgn{3})))
v   =   axis;
axis([min(t(2:end)*(10^time{3})),max(t(2:end)*(10^time{3})),v(3),v(4)])
figure(4)
plot(y(2:end)*(10^yn{3}),vg)
xlabel(sprintf('Centre position, (%sm)',yn{2}))
ylabel('Energy velocity, (c)')
