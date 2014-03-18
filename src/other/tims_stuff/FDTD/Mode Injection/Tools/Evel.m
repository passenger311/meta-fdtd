% Read file names
close all
clear all
clc
%tstart = tic;
% Read params

filename    =   'sim_params.txt';
fileID      =   fopen(filename,'r');
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
typ         =   in{1}{1};
typ         =   isequal(typ,'dx');
if typ == 1
    resolution  =   str2double(in{1}{2});
else
    res     =   str2double(in{1}{2});
end
foldername = 'Results';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end
tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xstart      =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xfinish     =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
central_lambda  =   str2double(in{1}{2});
fclose(fileID);


lightSI = 299792458;
%central_omega       =   1.4e15;
%central_freq        =   central_omega/(2*pi);
central_lambda      =   central_lambda*1e-6;

if  typ == 1
    resolution  =   resolution*1e-6;
else
    
    resolution          = central_lambda/res/(sqrt(11.68));
end
Courant             = 0.7;
TimeStepREAL        = (Courant*resolution)/lightSI;
% time                = 0;
% time                = round(time/TimeStepREAL);
% k0 = 2*pi/central_lambda;


filename    =   'evel.0.gpl';

fileID  =   fopen(filename,'r');
for loop = 1:8
    tline = fgetl(fileID);
end
[A,COUNTA] = textscan(fileID,'%n',1e9,'Delimiter','\n','CommentStyle','#');
fclose(fileID);
A = A{1};
n_time_steps = textscan(tline,'%s');
n_time_steps    =   str2double(n_time_steps{1}{end});
u.x  =   A(1:6:end);
u.y  =   A(2:6:end);
u.z  =   A(3:6:end);
o.x  =   A(4:6:end);
o.y  =   A(5:6:end);
o.z  =   A(6:6:end);
Unames  =   fieldnames(u);
x   =   1:length(u.x);


% to be removed later 
% filename = 'Angle.txt';
% fileID = fopen(filename);
% angle = str2double(fgetl(fileID));
% v_test = angle-40;
% fclose(fileID);
%remove


filename    =   'Results/evel_out.txt';

fileID  =   fopen(filename,'w');
for loopindex = 1:numel(Unames)
    y   =   u.(Unames{loopindex})';
    y2  =   o.(Unames{loopindex})';
    pos{loopindex} = y;
    disp{loopindex}     =   y2;
    fin = floor(length(x)/2);
    fit_params  =   least_squares(x(1:fin),y(1:fin));
    best.m  =   fit_params(1);
    error.m =   fit_params(2);
    best.c  =   fit_params(3);
    error.c =   fit_params(4);
    fit{loopindex}  =   (best.m.*x+best.c);
    fit_params2  =   least_squares(x(1:fin),y2(1:fin));
    best.m2  =   fit_params2(1);
    error.m2 =   fit_params2(2);
    best.c2  =   fit_params2(3);
    error.c2 =   fit_params2(4);
    fit2{loopindex}  =   (best.m2.*x+best.c2);
    f = figure('visible','off');
    plot(x*TimeStepREAL*n_time_steps*1e15,y*resolution*1e6,'b')
    hold on
    plot(x*TimeStepREAL*n_time_steps*1e15,fit{loopindex}*resolution*1e6,'r')
    hold off
    title(sprintf('Movement of centre of energy (%s)',Unames{loopindex}));
    xlabel('Time (fs)')
    ylabel('Centre position (um)');
    set(gcf,'color','white');
    legend('data','fit');
    filename    =   sprintf('Results/Centre_%s.eps',Unames{loopindex});
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(x*TimeStepREAL*n_time_steps*1e15,y2*resolution*1e6*2*sqrt(2*log(2)),'b')
    hold on
    plot(x*TimeStepREAL*n_time_steps*1e15,fit2{loopindex}*resolution*1e6*2*sqrt(2*log(2)),'r')
    hold off
    title(sprintf('Change in FWHM of pulse (%s)',Unames{loopindex}));
    xlabel('Time (fs)')
    ylabel('FWHM (um)');
    set(gcf,'color','white');
    legend('data','fit');
    filename    =   sprintf('Results/STD_%s.eps',Unames{loopindex});
    print(f,'-depsc','-r200','-painters',filename);
    vg.(Unames{loopindex})  =   best.m/n_time_steps/Courant;
    vg.er.(Unames{loopindex})   =   error.m/n_time_steps/Courant;
    d.(Unames{loopindex})   =   lightSI*best.m2/n_time_steps/Courant;
    d.er.(Unames{loopindex})   =   lightSI*error.m2/n_time_steps/Courant;
    
    % also remove later
    %vg.(Unames{1})  =   v_test;
    % REMOVE
    
    fprintf(fileID,'Energy velocity in %s direction = %g +/- %g\n',Unames{loopindex},vg.(Unames{loopindex}),vg.er.(Unames{loopindex}));
end
fprintf(fileID,'\n');
for loopindex = 1:numel(Unames)
    fprintf(fileID,'Energy dispersion in %s direction = %g +/- %g ms-1\n',Unames{loopindex},d.(Unames{loopindex}),d.er.(Unames{loopindex}));
end
fclose(fileID);

filename    =   'Results/Evel_fit.txt';

fileID  =   fopen(filename,'w');
time    =   x*TimeStepREAL*n_time_steps;
for loop = 1:length(time)
    fprintf(fileID,'%g\t%g\t%g\n',time(loop),pos{1}(loop),fit{1}(loop));
end
fclose(fileID);
filename    =   'Results/Evel_disp_fit.txt';

fileID  =   fopen(filename,'w');
time    =   x*TimeStepREAL*n_time_steps;
for loop = 1:length(time)
    fprintf(fileID,'%g\t%g\t%g\n',time(loop),disp{1}(loop),fit2{1}(loop));
end
fclose(fileID);

filename    =   'Results/Job_name.txt';
fileID  =   fopen(filename,'w');
a = cd;
ind = strfind(a,'/');
job_name = a(ind(end)+1:end);
fprintf(fileID,job_name);
fclose(fileID);

quit;





