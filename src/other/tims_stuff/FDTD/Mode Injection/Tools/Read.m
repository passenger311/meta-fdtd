% Read file names
close all
clear all
clc
tstart = tic;
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
time                = 0;
time                = round(time/TimeStepREAL);
k0 = 2*pi/central_lambda;
% cd 28388.localhost/Data/
list = dir('Hz_*.gpl');
num = max(size(list));

for n = 1:num
    name = list(n).name;
    ind1 = strfind(name,'_');
    ind2 = strfind(name,'.0');
    str_num = str2num(name(ind1+1:ind2-1));
    num_list(n,1) = {str_num};

    num_list(n,2) = {(name(ind1+1:ind2-1))};
end

num_list = sortrows(num_list,1);
ind = [num_list{:,1}] > time;
num_list = num_list(ind,:);

num = max(size(num_list));
per = 0;
fprintf('Fitting ... %0.1f%%\n',per);

for n = 1:num
    clc
    pern     =   floor(1000*n/num);
    pern     =   pern/10;
    
    if pern>per
        
        fprintf('Fitting ... %0.1f%%\n',pern);
        per = pern;
        
    end
    
    ncyc = num_list{n,2};
    filename = sprintf('Hz_%s.0.gpl',ncyc);
    
    file = fopen(filename,'r');
    [A,COUNTA] = textscan(file,'%n',1e5,'Delimiter','\n','CommentStyle','#');
    fclose(file);
    A = A{1};
    
    ind  = isnan(A) == 0;
    A = A(ind);
    A_full = A;
    A = A(xstart:xfinish);
    C = xstart:1:xfinish;
    if n == 1
        a = (abs(min(A))+abs(max(A)))/2;
        A_s = sign(A);
        A_s2 = A_s(1:end-1)+A_s(2:end);
        ind = find(A_s2 == 0);
        ind2 = ind(2:end)-ind(1:end-1);
        ind2 = ind2(4:end-2);
        l_2 = mean(ind2);
        lambda = l_2*2;
        freq = 1/lambda;
        
        w = 2*pi*freq;
        phi = 0;
        b = (find(A == max(A)) + find(A == min(A)))/2;
        b = C(round(mean(b)));
        ind = find(A >= a/2);
        fwhm = ind(end)-ind(1);
        c = fwhm/(2*sqrt(2*log(2)));
        p(1) = a;
        p(2) = w;
        p(3) = phi;
        p(4) = b;
        p(5) = c;
%          p(6) = a/2;
%          p(7) = w+0.02;
%          p(8) = phi;
%          p(9) = b-100;
%          p(10) = c;
    else
        p(1) = amp_mat(n-1);
        p(2) = omega_mat(n-1);
        p(3) = phase_mat(n-1);
        p(4) = centre_mat(n-1);
        p(5) = width_mat(n-1);
    end
% %      cd ..
% %      cd ..
    s_params = p;
    options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
    estimate = fminsearch(@myfit,s_params,options,C',A);
    
    A_mat(n,:) = A;
    amp_mat(n) = estimate(1);
    omega_mat(n) = estimate(2);
    phase_mat(n) = estimate(3);
    centre_mat(n) = estimate(4);
    width_mat(n) = estimate(5);
%     amp_mat2(n) = estimate(6);
%     omega_mat2(n) = estimate(7);
%     phase_mat2(n) = estimate(8);
%     centre_mat2(n) = estimate(9);
%     width_mat2(n) = estimate(10);
    envel_mat(n,:) = estimate(1).*exp(-(((C-estimate(4)).^2)./(2*(estimate(5)^2))));
%     if n > 1
%         
%         diff(1) = abs((amp_mat(n)/amp_mat(n-1))-1);
%         diff(2) = abs((omega_mat(n)/omega_mat(n-1))-1);
%         diff(3) = abs((centre_mat(n)/centre_mat(n-1))-1);
%         diff(4) = abs((width_mat(n)/width_mat(n-1))-1);
%         
%         m_diff = max(diff);
%         
%         if m_diff > 0.5
%             
%             figure(1)
%             plot(A_full);
%             hold on
%             D = amp_mat(end).*sin((omega_mat(end).*C)+phase_mat(end)).*exp(-((C-centre_mat(end)).^2)./(2*(width_mat(end)^2)));
%             plot(C',D,'r');
%             %   axis([C(1) C(end) -0.8 0.8])
%             hold off
%             
%             pause(10)
%         end
%     end
%     cd 28388.localhost/Data/
end
%      cd ..
%      cd ..
foldername = 'Results';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end
%  x = A;
%  N=length(x);
%  Fs = 1;
% %this part of the code generates that frequency axis
% if mod(N,2)==0
%     k=-N/2:N/2-1; % N even
% else
%     k=-(N-1)/2:(N-1)/2; % N odd
% end
% T=N/Fs;
% freq=k/T;  %the frequency axis
% % D = amp_mat(end).*sin((omega_mat(end).*C)+phase_mat(end)).*exp(-((C-centre_mat(end)).^2)./(2*(width_mat(end)^2)));
% %takes the fft of the signal, and adjusts the amplitude accordingly
% X=fft(x)/N; % normalize the data
% X=fftshift(X); %shifts the fft data so that it is centered
%  y = x;
% [YfreqDomain,freq] = centeredFFT(y,Fs);
% figure(2)
% plot(freq*2*pi,abs(YfreqDomain))
% zeroPadFactor = nextpow2(length(y)) + 3;
% [a,b] = positiveFFT_zero_padding(y,Fs,2^zeroPadFactor);
% figure(3)
% b   =   (2*pi)./(resolution./b);
% plot(b,abs(a))
% aba = abs(a);
% mx = max(aba);
% ind = find(aba == mx);
% f = b(ind);
% lam_h = 1/f;
% lam_h = lam_h*resolution;
% k_h = 2*pi/lam_h;
% neff_h = k_h/k0;
% ind = find(aba >= mx/2);
% fwhm = b(ind(end))-b(ind(1));
% c = fwhm/(2*sqrt(2*log(2)));
% p = [];
% p(1) = mx;
% p(2) = f;
% p(3) = c;
% s_params = p;
% options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
% estimate = fminsearch(@myfit2,s_params,options,b',aba);
% f   =   estimate(2);
% lam_h = 1/f;
% lam_h = lam_h*resolution;
% k_h = 2*pi/lam_h;
% neff_h = k_h/k0;
%     fitted_curve = estimate(1).*exp(-((b-estimate(2)).^2)./(2*(estimate(3)^2)));
%     figure(4)
%     plot(b*2*pi,abs(a))
%     hold on
%     plot(b*2*pi,fitted_curve,'r')

%f = figure('visible','off');
%figure(1)
%plot(C',A);
%hold on
  D = (amp_mat(end)).*sin((omega_mat(end).*C)+phase_mat(end)).*exp(-((C-centre_mat(end)).^2)./(2*(width_mat(end)^2)));
% %D = D+(amp_mat(end)/2).*sin(((omega_mat(end)).*C)+phase_mat(end)).*exp(-((C-(centre_mat(end)+20)).^2)./(2*(width_mat(end)^2)));
% %D = D+amp_mat2(end).*sin((omega_mat2(end).*C)+phase_mat2(end)).*exp(-((C-centre_mat2(end)).^2)./(2*(width_mat2(end)^2)));
% %D = D + amp_mat(end-1)/10.*sin((omega_mat(end-1).*C)+phase_mat(end-1)).*exp(-((C-centre_mat(end-1)).^2)./(2*(width_mat(end-1)^2)));
% plot(C',D,'r');
% hold on
% E = amp_mat(end).*exp(-(((C-centre_mat(end)).^2)./(2*(width_mat(end)^2))));
% plot(C',E,'k');
% hold on
% plot(centre_mat(end)*ones(1,2),[-1 1],'k')
% hold on
% plot(centre_mat(end)*ones(1,2)+width_mat(end),[-1 1],'k')
% hold on
% plot(centre_mat(end)*ones(1,2)-width_mat(end),[-1 1],'k')
%hold off
% title('Fit comparison')
% xlabel('X position (grid cells)');
% ylabel('Hz field (arb)');
%axis([C(1) C(end) -max(abs(A))*1.1 max(abs(A))*1.1]);
% filename = 'Results/Fit_comparison.eps';
%print(f,'-depsc','-r200','-painters',filename);
%return
% figure(2)
% plot(C',(A-D'))

num_list = [num_list{:,1}];
num_list = num_list';
time = num_list*TimeStepREAL;
%f = figure('visible','off');
% figure(2)
% plot(time,amp_mat,'Linewidth',2);
% ax = axis;
% axis([time(1) time(end) 0 ax(end)]);
% title('Amplitude decay');
% xlabel('Time (s)');
% ylabel('Hz field amplitude (arb)');
% filename = 'Results/Amplitude.eps';
%print(f,'-depsc','-r200','-painters',filename);
%f = figure('visible','off');
% figure(3)
% plot(time,omega_mat,'Linewidth',2);
% ax = axis;
% axis([time(1) time(end) ax(3) ax(4)]);
% title('Wavevector evolution');
% xlabel('Time (s)');
% ylabel('Width (grid cells)');
% filename = 'Results/Wavevector.eps';
%print(f,'-depsc','-r200','-painters',filename);

%f = figure('visible','off');
% figure(4)
% plot(time,centre_mat,'Linewidth',2);
% ax = axis;
% axis([time(1) time(end) ax(3) ax(4)]);
% title('Centre position movement');
% xlabel('Time (s)');
% ylabel('Centre position (grid cells)');
% filename = 'Results/Centre_position.eps';
%print(f,'-depsc','-r200','-painters',filename);

%f = figure('visible','off');
% figure(5)
% plot(time,width_mat,'Linewidth',2);
% ax = axis;
% axis([time(1) time(end) ax(3) ax(4)]);
% title('Width evolution');
% xlabel('Time (s)');
% ylabel('Width (grid cells)');
% filename = 'Results/Width.eps';
%print(f,'-depsc','-r200','-painters',filename);

%f = figure('visible','off');
% figure(6)
% plot(time,log(amp_mat),'Linewidth',2);
% ax = axis;
% axis([time(1) time(end) ax(3) ax(4)]);
% title('Log amplitude');
% xlabel('Time (s)');
% ylabel('Log(amplitude) (arb)');
% filename = 'Results/Log(amp).eps';
%print(f,'-depsc','-r200','-painters',filename);
k0 = 2*pi/central_lambda;




fprintf('\n\nSaving results to txt files.\n');


outfile = 'Results/fit_params.txt';
fileID = fopen(outfile,'w');
fprintf(fileID,'Time(s)\tAmplitude\tOmega\tPhase\tCentre\tWidth\n');

for n = 1:length(time)
    fprintf(fileID,'%g\t%g\t%g\t%g\t%g\t%g\n',time(n),amp_mat(n),omega_mat(n),phase_mat(n),centre_mat(n),width_mat(n));
end

fclose('all');

out_file = 'Results/output.txt';
fileID = fopen(out_file,'w');
fprintf(fileID,'Time(s)\tAmplitude\tNeff\tCentre\tFWHM\tlog(amp)\n');
for n = 1:length(time)
    lam_l = 2*pi/omega_mat(n);
    lam_l = lam_l*resolution;
    k_l = 2*pi/lam_l;
    Neff = k_l/k0;
    fprintf(fileID,'%g\t%g\t%g\t%g\t%g\t%g\n',time(n),amp_mat(n),Neff,centre_mat(n)*resolution,(2*sqrt(2*log(2)))*width_mat(n)*resolution,log(amp_mat(n)));
end

fclose('all');

out_file = 'Results/fit.txt';
fileID = fopen(out_file,'w');
fprintf(fileID,'x(dx)\tReal_data\tFit_data\n');
for n = 1:length(C)
    x = C(n);
    real_y = A(n);
    fit_y = D(n);
    
    fprintf(fileID,'%g\t%g\t%g\n',x,real_y,fit_y);
end

fclose('all');


fprintf('\nCalculating derived parameters\n');
xstart = 480;%input('\n\nxstart = ');

centre = centre_mat(xstart:end);
log_amp = log(amp_mat(xstart:end));
time = time(xstart:end);
FWHM_mat = (2*sqrt(2*log(2))).*width_mat(xstart:end).*resolution;



% Vg calculation
m_centre = mean((centre(2:end)-centre(1:end-1))./TimeStepREAL);
c_centre = mean(centre'-(m_centre.*time));
s_params = [m_centre, c_centre];
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
estimate = fminsearch(@linfit,s_params,options,time,centre');
m_c   =   estimate(1);
c_c   =   estimate(2);



%f = figure('visible','off');
% figure(7)
% plot(time,centre);
% hold on
centre_fit = (estimate(1).*time)+estimate(2);

centre_diff = centre-centre_fit';
c_new   =   max(abs(centre_diff))/2;
ex = 0;
loop = 0;
while ex == 0
    loop = loop +1;
    ind =   find(centre_diff>=-c_new & centre_diff<=c_new);
    nm  =   length(ind)/length(centre);
    if (nm >= 0.68 & nm <=0.69);
        ex = 1;
    elseif nm < 0.68
        c_new = c_new*1.001;
    elseif nm > 0.69
        c_new = c_new*0.999;
    end
    if loop == 5000
        break
    end
end
        
par_u   =   centre_fit+c_new;
par_l   =   centre_fit-c_new;
dx  =   time(end)-time(1);
dy  =   par_u(end)-par_l(1);
m_u =   dy/dx;
c_u =   par_l(1)-time(1)*m_u;
dy  =   par_l(end)-par_u(1);
m_l =   dy/dx;
c_l =   par_u(1)-time(1)*m_l;
fit_l   =   (m_l*time) + c_l;
fit_u   =   (m_u*time)+c_u;
% 
% plot(time,centre_fit,'r')
% hold on
% plot(time,par_u,'k')
% hold on
% plot(time,par_l,'k')
% hold on
% plot(time,fit_l,'g')
% hold on
% plot(time,fit_u,'g')
% hold off
% title('Vg calculation');
% xlabel('Time (s)')
% ylabel('Centre position (m)');
% filename = 'Results/Vg_calc.eps';
%print(f,'-depsc','-r200','-painters',filename);
vg = (estimate(1)*resolution)/lightSI;
vg_u = (m_u*resolution)/lightSI;
vg_l = (m_l*resolution)/lightSI;
vg_er   =   abs(vg-vg_u);





% Dispersion calculation
m_disp = mean((FWHM_mat(2:end)-FWHM_mat(1:end-1))./TimeStepREAL);
c_disp = mean(FWHM_mat'-(m_disp.*time));
s_params = [m_disp, c_disp];
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
estimate = fminsearch(@linfit,s_params,options,time,FWHM_mat');
%f = figure('visible','off');
% figure(8)
% plot(time,FWHM_mat);
% hold on
disp_fit = (estimate(1).*time)+estimate(2);


disp_diff = FWHM_mat-disp_fit';
d_new   =   max(abs(disp_diff))/2;
ex = 0;
loop = 0;
while ex == 0
    loop = loop +1;
    ind =   find(disp_diff>=-d_new & disp_diff<=d_new);
    nm  =   length(ind)/length(FWHM_mat);
    if (nm >= 0.68 & nm <=0.69);
        ex = 1;
    elseif nm < 0.68
        d_new = d_new*1.001;
    elseif nm > 0.69
        d_new = d_new*0.999;
    end
    if loop == 5000
        break
    end
    
end
        
par_u   =   disp_fit+d_new;
par_l   =   disp_fit-d_new;
dx  =   time(end)-time(1);
dy  =   par_u(end)-par_l(1);
m_u =   dy/dx;
c_u =   par_l(1)-time(1)*m_u;
dy  =   par_l(end)-par_u(1);
m_l =   dy/dx;
c_l =   par_u(1)-time(1)*m_l;
fit_l   =   (m_l*time) + c_l;
fit_u   =   (m_u*time)+c_u;


% plot(time,disp_fit,'r')
% hold on
% plot(time,par_u,'k')
% hold on
% plot(time,par_l,'k')
% hold on
% plot(time,fit_l,'g')
% hold on
% plot(time,fit_u,'g')
% hold off
% title('Disp calculation');
% xlabel('Time (s)')
% ylabel('FWHM (m)');
% filename = 'Results/Disp_calc.eps';
%print(f,'-depsc','-r200','-painters',filename);
disp = estimate(1);
disp_u = m_u;
disp_l = m_l;
disp_er   =  abs(disp-disp_u);


% Loss calculation
m_amp = mean((log_amp(2:end)-log_amp(1:end-1))./TimeStepREAL);
c_amp = mean(log_amp'-(m_amp.*time));
s_params = [m_amp, c_amp];
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
estimate = fminsearch(@linfit,s_params,options,time,log_amp');
%f = figure('visible','off');
% figure(9)
% plot(time,log_amp);
% hold on
amp_fit = (estimate(1).*time)+estimate(2);

amp_diff = log_amp-amp_fit';
a_new   =   max(abs(amp_diff))/2;
ex = 0;
loop = 0;
while ex == 0
    loop = loop +1;
    ind =   find(amp_diff>=-a_new & amp_diff<=a_new);
    nm  =   length(ind)/length(log_amp);
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
        
par_u   =   amp_fit+a_new;
par_l   =   amp_fit-a_new;
dx  =   time(end)-time(1);
dy  =   par_u(end)-par_l(1);
m_u =   dy/dx;
c_u =   par_l(1)-time(1)*m_u;
dy  =   par_l(end)-par_u(1);
m_l =   dy/dx;
c_l =   par_u(1)-time(1)*m_l;
fit_l   =   (m_l*time) + c_l;
fit_u   =   (m_u*time)+c_u;



% plot(time,amp_fit,'r')
% hold on
% plot(time,par_u,'k')
% hold on
% plot(time,par_l,'k')
% hold on
% plot(time,fit_l,'g')
% hold on
% plot(time,fit_u,'g')
% hold off
% title('Loss calculation');
% xlabel('Time (s)')
% ylabel('Log(amp) (m)');
% filename = 'Results/Loss_calc.eps';
%print(f,'-depsc','-r200','-painters',filename)
alpha = estimate(1)*2;
alpha_u = m_u*2;
alpha_l = m_l*2;
alpha_er   =   abs(alpha-alpha_u);




lam_l = 2*pi/omega_mat(xstart);
lam_l = lam_l*resolution;
k_l = 2*pi/lam_l;

neff_l = k_l/k0;

lam_h = 2*pi/omega_mat(end);
lam_h = lam_h*resolution;
k_h = 2*pi/lam_h;
neff_h = k_h/k0;
fprintf('\nSaving parameters.\n');
filename = 'Results/parameters.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'Central Wavelength = %g m\n',central_lambda);
fprintf(fileID,'dx = %g m\n',resolution);
fprintf(fileID,'dt = %g s\n',TimeStepREAL);
fprintf(fileID,'Neff = %g to %g \n',neff_l,neff_h);
fprintf(fileID,'Group Velocity = %g +/- %g (c) \n',vg,vg_er);
fprintf(fileID,'Loss = %g +/- %g (1/s) \n',alpha,alpha_er);
fprintf(fileID,'Dispersion = %g +/-%g (m/s)',disp,disp_er);

fclose('all');
%f = figure('visible','off');
% figure(10)
% plot(centre_mat,amp_mat,'LineWidth',2);
% title('Spatial decay');
% xlabel('Centre_position');
% ylabel('Amplitude');
% filename = 'Results/Spatial_decay.eps';
%print(f,'-depsc','-r200','-painters',filename);
% ind = 0;
% inc = floor(num/3);
% col = 'brgkcy';
% in = 0;
% for ind = 1:inc-1:num
%     in = in+1;
%     hold on
%     plot(C',A_mat(ind,:),col(in))
%     hold on
%     plot(C',envel_mat(ind,:),col(in))
% end
% hold off


filename = 'Results/Pulse_decay.txt';
fileID = fopen(filename,'w');
ind1 = max(amp_mat)/2;
ind1 = abs(amp_mat-ind1);
ind1 = find(ind1 == min(ind1));
ind2 = max(amp_mat)/4;
ind2 = abs(amp_mat-ind2);
ind2 = find(ind2 == min(ind2));
ind3 = max(amp_mat)/10;
ind3 = abs(amp_mat-ind3);
ind3 = find(ind3 == min(ind3));

for loop = 1:length(C)
    x = C(loop);
    p1  =   A_mat(1,loop)./amp_mat(1);
    e1  =   envel_mat(1,loop)./amp_mat(1);
    p2  =   A_mat(ind1,loop)./amp_mat(1);
    e2  =   envel_mat(ind1,loop)./amp_mat(1);
    p3  =   A_mat(ind2,loop)./amp_mat(1);
    e3  =   envel_mat(ind2,loop)./amp_mat(1);
    p4  =   A_mat(ind3,loop)./amp_mat(1);
    e4  =   envel_mat(ind3,loop)./amp_mat(1);
    fprintf(fileID,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',x,p1,e1,p2,e2,p3,e3,p4,e4);
end


fclose('all');

filename = 'Results/Pulse_decay2.txt';
fileID = fopen(filename,'w');

for loop = 1:length(amp_mat)
    x = centre_mat(loop);
    amp =   amp_mat(loop)/amp_mat(1);
    fprintf(fileID,'%g\t%g\n',x,amp);
end

fprintf('\nDone\n');
t   =   toc(tstart);
fprintf('\nTime taken = %0.1fs\n',t);
fclose('all');
quit;
