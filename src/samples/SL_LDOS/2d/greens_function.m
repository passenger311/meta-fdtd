%%% Clear
clear all;
clf

%%% Setup parameters
format long;
resolution=10;
DT=0.7;
ncyc_max=8*(32767+1);
nstep=1;
natt=9000+1;
nhwhm=330;
l_inv=0.0066666666666667;

%%% Set up current source
ncyc=(1:nstep:ncyc_max)';
ncyc_full=(1:ncyc_max)';
om=2*pi*l_inv;
g=sqrt(log(2))/(nhwhm*DT);

j=cos(om.*(ncyc-natt)*DT).*exp(-g.^2.*((ncyc-natt)*DT).^2);
j_full=cos(om.*(ncyc_full-natt)*DT).*exp(-g.^2.*((ncyc_full-natt)*DT).^2);
%plot(ncyc,j);break;

%%% Print current source to file
%fid1=fopen('j.dat','w');
%fprintf(fid1,'  %e\n', j_full);
%fclose(fid1);

%%% Calculate polarisation
p_full=j_full;
for k=2:size(p_full,1);
  p_full(k)=p_full(k-1)+j_full(k);
end;
p_full=p_full*DT;
%p2_full=cumtrapz(j_full)*DT; %%% Produces wrong results as meta uses above integration scheme
p=p_full(ncyc);
%p2=p2_full(ncyc);
%plot(ncyc,p,ncyc,p2);break;


%%% Perform Fourier Transform
p_om=fft(p);
%p2_om=fft(p2);
%j_om=fft(j);

%%% Normalise
%p_om=p_om/(size(p_om,1)/2)/pi;
%p2_om=p2_om/(size(p2_om,1)/2)/pi;
%j_om=j_om/(size(j_om,1)/2)/pi;

%%% Calculate frequency range
freq=1./DT*(ncyc-1)/(ncyc_max)/nstep;

%%% Get rid of high frequency part (which is complex conjugate of low freqency part)
freq(size(freq,1)/2+2:size(freq,1)) = [];
p_om(size(p_om,1)/2+2:size(p_om,1)) = [];
%p2_om(size(p2_om,1)/2+2:size(p2_om,1)) = [];
%j_om(size(j_om,1)/2+2:size(j_om,1)) = [];

%%% Get rid of zero frequency
freq(1) = [];
p_om(1) = [];
%p2_om(1) = [];
%j_om(1) = [];

%%% Load signal electric field from meta
filename=sprintf('E.0.gpl');
fprintf(1,'\nloading file: %s\n', filename);
fid1 = fopen(filename, 'r');
if (fid1==-1) continue; end;
for k=1:13
  tline = fgetl(fid1);
end
tmp = textscan(fid1, '%n %n %n');
fclose(fid1);
tmp2=zeros(size(tmp{1},1),size(tmp,2));
for k=1:3
  tmp2(:,k)=tmp{k};
end
e=tmp2(1:size(p),:);
clear tmp tmp2;

%%% Perform Fourier Transform, normalise and clip zero and high frequencies
e_om=fft(e);
%e_om=e_om/(size(e_om,1)/2)/pi;
e_om(size(e_om,1)/2+2:size(e_om,1),:) = [];
e_om(1,:) = [];

%%% Some tests
%[a,b]=max(abs(p_om));
%c=freq(b);
%[a1,b1]=max(abs(ex_om));
%c1=freq(b1);
%fprintf('amplitude difference: %e\n',a/a1)
%fprintf('frequencies: %e %e\t at points %i %i\n',c,c1,b,b1)
%fprintf('frequency diff: %e %e\n',freq(2)-freq(1))

%%% Reduce size of vectors
%plot(abs(p_om))
[tmp, min_f]=min(abs(abs(p_om)-max(abs(p_om))./10));
[tmp, pmax_f]=max(abs(p_om));
if min_f<pmax_f
  [tmp, max_f]=min(abs(abs(p_om(min_f+10:size(p_om)))-max(abs(p_om))./10));
  max_f=min_f+10+max_f;
else
  max_f=min_f;
  [tmp, min_f]=min(abs(abs(p_om(1:max_f-10))-max(abs(p_om))./10));
end
p_om([1:min_f max_f:size(p_om,1)]) = [];
e_om([1:min_f max_f:size(e_om,1)],:) = [];
freq([1:min_f max_f:size(freq,1)]) = [];
%p2_om([1:min_f max_f:size(p2_om,1)]) = [];
%j_om([1:min_f max_f:size(j_om,1)]) = [];


%%% Plotting j and p2 for checks (have to be uncommented above!)
%plot(freq,j_om,freq,i.*2.*pi.*freq.*p2_om,'x',freq,i.*2.*pi.*freq.*p_om,'o');break;
%plot(freq,real(i.*2.*pi.*freq.*p_om./j_om),freq,real(i.*2.*pi.*freq.*p2_om./j_om),'o');break;
%plot(freq,angle(p_om),freq,angle(p2_om),'x',freq,angle(j_om./(i.*2.*pi.*freq)),'o');break;
%plot(freq,angle(p2_om)-angle(p_om),freq,angle(j_om./(i.*2.*pi.*freq))-angle(p_om),'o');break;

%%% Calculate Green's function
g=zeros(size(e_om));
for k = 1:3
  g(:,k)=imag(e_om(:,k)./p_om);
end;
%g2=imag(ex_om./p2_om);
%g3=imag(i*2*pi*freq.*ex_om./j_om);

%%% Change to SI units
freq_real=freq*2.9979e5/resolution;
g_real=g*(1e9/resolution)^2;
%g2_real=g2*(1e9/resolution)^2;
%g3_real=g3*(1e9/resolution)^2;
g_theo=pi^2./2.*freq.^2*(1e9/resolution)^2;%4.*pi^2./3.*freq.^2*(1e9/resolution)^2;

figure(1);
subplot(2,2,1);
plot(freq_real, g_real(:,1),'o', freq_real, g_theo, 'red', 'LineWidth', 1.5)
xlabel('Frequency [THz]','FontSize',15);
ylabel('Greens function [1/m^2]','FontSize',15);
set(gca,'FontSize',15);
%print(gcf, '-dpng', '-r150', 'greens_function_x.png');

subplot(2,2,2);
plot(freq_real, g_real(:,2),'o', freq_real, g_theo, 'red', 'LineWidth', 1.5)
xlabel('Frequency [THz]','FontSize',15);
ylabel('Greens function [1/m^2]','FontSize',15);
set(gca,'FontSize',15);
%print(gcf, '-dpng', '-r150', 'greens_function_y.png');

subplot(2,2,3);
plot(freq_real, g_real(:,3),'o', freq_real, g_theo, 'red', 'LineWidth', 1.5)
xlabel('Frequency [THz]','FontSize',15);
ylabel('Greens function [1/m^2]','FontSize',15);
set(gca,'FontSize',15);
%print(gcf, '-dpng', '-r150', 'greens_function_z.png');
print(gcf, '-dpng', '-r150', 'greens_function.png');

%%% Fit function
%F = @(x,xdata)x(1)+x(2).*(xdata-x(3)).^3;
%x0 = [ 0 5e11 0 ];
%[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,freq_real,g_real);
%fprintf('Fit value: %e\n', x);
%plot(freq_real,F(x,freq_real),'green')
%hold off;

%%% Print Green's function to file
fid1=fopen('greens_function.dat','w');
fprintf(fid1,'%e %e %e %e %e\n', [freq_real g_real g_theo]');
fclose(fid1);

figure(2)
plot(freq_real, g_real(:,1)./g_theo, freq_real, g_real(:,2)./g_theo, freq_real, g_real(:,3)./g_theo, 'LineWidth', 1.5)
xlabel('Frequency [THz]','FontSize',15);
ylabel('Enhancement','FontSize',15);
set(gca,'FontSize',15);
print(gcf, '-dpng', '-r150', 'enhancement.png');

fid1=fopen('enhancement.dat','w');
fprintf(fid1,'%e %e %e %e\n', [freq_real g_real(:,1)./g_theo g_real(:,2)./g_theo g_real(:,3)./g_theo]');
fclose(fid1);

exit
