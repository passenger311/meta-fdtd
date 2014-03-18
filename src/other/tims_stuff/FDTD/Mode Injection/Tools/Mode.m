function Mode(dx,dt)

% Constants
si.c    =   299792458;      % Speed of light in vacuum (m/s)
dt_r    =   dt*dx*1e-6/si.c;     % Time step size

% Create results folder

foldername = 'Results';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end


% Step 1: Plot the mode profiles

filename    =   'modeplot_stability_1.set';
fileID      =   fopen(filename,'r');

% Read array dimension from file

for a = 1:3
    tline = fgetl(fileID);
    if a == 2
        dim = str2num(tline);
    end
end

xst     =   dim(1);
xfi     =   dim(2);
inc     =   dim(3);
xle     =   ceil((xfi-xst)/inc)+1;
yst     =   dim(4);
yfi     =   dim(5);
inc     =   dim(6);
yle     =   ceil((yfi - yst)/inc);

% Read out array

fields  =   zeros(xle*yle,12);      % Initalise array
loop    =   0;

while 1
    loop = loop+1;
    if isempty(tline) == 0
        fields(loop,:) = str2num(tline);
    end
    tline = fgetl(fileID);
    if ~ischar(tline),   break, end
    if isequal(tline,')SET'), break,end
end

fclose(fileID);

% Find size of gain region from .in file

filename    =   'geo_geo_gain.in';
fileID      =   fopen(filename,'r');

for a = 1:4
    tline = fgetl(fileID);
end
fclose(fileID);

tline   =   str2num(tline);
gst     =   tline(1);
gfi     =   tline(2);
a       =   1:yle;
b       =   1:xle;
a       =   a*2*dx*1e6;
b       =   b*2*dx*1e6;
nam     =   {'Ex','Ey','Ez','Hx','Hy','Hz'};

% Create Mode_profile folder

foldername = 'Results/Mode_profiles';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end

cl = 'k';
wid = 1;
ind = [1 3 11];

for loop = 1:3
    n = nam{ceil(ind(loop)/2)};
    f1 = fields(:,ind(loop));
    f1 = reshape(f1,xle,yle);
    f1 = f1';
    f = figure('visible','off');
    pcolor(b,a,f1)
    shading interp
    hold on
    x = [0 xle*2]*dx*1e6;
    y = [62 62]*dx*1e6;
    plot(x,y,cl,'LineWidth',wid)
    hold on
    y = [92 92]*dx*1e6;
    plot(x,y,cl,'LineWidth',wid)
    hold on
    y = [141 141]*dx*1e6;
    plot(x,y,cl,'LineWidth',wid)
    hold on
    y = [62 92]*dx*1e6;
    x = [gst gst]*dx*1e6;
    plot(x,y,cl,'LineWidth',wid)
    hold on
    x = [gfi gfi]*dx*1e6;
    plot(x,y,cl,'LineWidth',wid)
    title(n,'FontSize',16)
    xlabel('X position (um)','FontSize',14)
    ylabel('Y position (um)','FontSize',14)
    colorbar
    filename = sprintf('%s/%s.eps',foldername,n);
    print(f,'-depsc','-r200','-painters',filename);
end


% Step 2: Plot the population of the 4 energy levels

filename    =   'dens_pump.0.gpl';
fileID      =   fopen(filename,'r');

loop = 0;

for a = 1:13
    tline = fgetl(fileID);
    if a == 8
        dim2 = str2num(tline(3:end));
    end
end

dens_mat    =   zeros(round(dim2(2)/dim2(3))+1,4);

while 1
    loop = loop+1;
    tline = fgetl(fileID);
    if ~ischar(tline),   break, end
    
    if isempty(tline) == 0
        dens_mat(loop,:) = str2num(tline);
    end
end
fclose(fileID);

time    =   1:length(dens_mat);
x_axis  =   time*dt_r*dim2(3);
x_axis  =   x_axis*1e12;
xst     =   0;
xfi     =   max(x_axis);
dens_mat(:,5) = dens_mat(:,3)-dens_mat(:,2);
a       =   round(length(dens_mat)/2);
av_E(1) =   mean(dens_mat(a:end,1));
av_E(2) =   mean(dens_mat(a:end,2));
av_E(3) =   mean(dens_mat(a:end,3));
av_E(4) =   mean(dens_mat(a:end,4));
av_E(5) =   mean(dens_mat(a:end,5));

y   =   [1 1];
ti  =   {'N0' 'N1' 'N2' 'N3' 'N2-N1'};

foldername = 'Results/Population Density';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end

for loop = 1:5
    f = figure('visible','off');
    plot(x_axis,dens_mat(:,loop))
    hold on
    plot([xst xfi],y*av_E(loop),'k')
    hold off
    axis([xst xfi 0 1])
    ab = sprintf('%g',av_E(loop));
    text(xfi*0.85,av_E(loop)+0.04,ab,'FontSize',12)
    title(ti{loop},'FontSize',16)
    xlabel('Time (ps)','FontSize',14)
    ylabel('Population density','FontSize',14)
    filename = sprintf('%s/%s.eps',foldername,ti{loop});
    print(f,'-depsc','-r200','-painters',filename);
end

% Step 3: Plot the time signal and frequency spectrum

filename    =   'Point_Hz.0.gpl';
fileID      =   fopen(filename,'r');


loop = 0;
for a = 1:13
    tline = fgetl(fileID);
    if a == 8
        dim2 = str2num(tline(3:end));
    end
end

time_sig    =   zeros(round(dim2(2)/dim2(3))+1,1);

while 1
    loop = loop+1;
    tline = fgetl(fileID);
    if ~ischar(tline),   break, end
    if isempty(tline) == 0
        time_sig(loop,:) = str2num(tline);
    end
end
fclose(fileID);


x   =   1:length(time_sig);
x   =   x*dt_r*dim2(3)*1e12;

f1  =   figure('visible','off');
plot(x,time_sig)
yax =   axis;
axis([0 max(x) yax(3:end)])
title('Hz field time signal','FontSize',16)
xlabel('Time (ps)','FontSize',14)
ylabel('Hz','FontSize',14)

sttime  =   6e-12;

time_sig    =   time_sig(round(sttime/(dt_r*dim2(3))):end);
x   =   1:length(time_sig);
y   =   time_sig;
x   =   x*dt_r*dim2(3)*1e12;

zeroPadFactor = nextpow2(length(time_sig)) + 3;

Fs      =   1/(dt_r*dim2(3));
L       =   length(time_sig);
NFFT    =   2^nextpow2(L);
Y       =   fft(y,NFFT)/L;
f       =   Fs/2*linspace(0,1,NFFT/2+1);
f       =   si.c./f;
amp     =   2*abs(Y(1:NFFT/2+1));

N       =   2^zeroPadFactor;
k       =   0:N-1;                          %create a vector from 0 to N-1
T       =   N/Fs;                           %get the frequency interval
freq    =   k/T;                            %create the frequency range
X       =   fft(y,N)/length(time_sig);      % normalize the data

%only want the first half of the FFT, since it is redundant
cutOff  =   ceil(N/2);

%take only the first half of the spectrum
X       =   X(1:cutOff);
freq    =   freq(1:cutOff);
f2      =   Fs/2*linspace(0,1,N/2+1);
f2      =   f2(1:end-1);
f       =   si.c./f2;
amp     =   2*abs(X(1:N/2));

% Find peaks
[pks,locs]  =   findpeaks(amp,'MINPEAKHEIGHT',max(amp)/10,'MINPEAKDISTANCE',100);


% Plot single-sided amplitude spectrum.
f2  =   figure('visible','off');
plot(f*1e6,amp)
hold on
plot(f(locs)*1e6,pks*1.005,'k^','markerfacecolor',[1 0 0]);
hold off
for loop = 1:length(locs)
    lab = sprintf('%g',f(locs(loop))*1e6);
    text(f(locs(loop))*1e6+0.014,pks(loop),lab,'FontSize',14);
end
title('Single-Sided Amplitude Spectrum of y(t)','FontSize',16)
xlabel('Wavelength (um)','FontSize',14)
ylabel('|Y(f)|','FontSize',14)
axis([1.2 1.7 0 1.1*max(amp)])
foldername = 'Results/Frequency Spectrum';
if (exist(foldername, 'dir') == 0)
    mkdir(foldername)
end

filename = sprintf('%s/Hz_Timesignal.eps',foldername);
print(f1,'-depsc','-r200','-painters',filename);
filename = sprintf('%s/FT_Hz_Timesignal.eps',foldername);
print(f2,'-depsc','-r200','-painters',filename);

% Fit gaussian and calculate Q factor
f2  =   figure('visible','off');
plot(f*1e6,amp)
hold on
for n = 1:2
    hw = 100;
    st = locs(n)-hw;
    fi = locs(n)+hw;
    tmp_f = f(st:fi)*1e6;
    tmp_amp = amp(st:fi);
    a = amp(locs(n));
    b = f(locs(n))*1e6;
    tmp_c = tmp_f(tmp_amp>=(a/2));
    FWHM = (max(tmp_c)-min(tmp_c));
    c = FWHM/(2*sqrt(2*log(2)));
    p(1) = a;
    p(2) = b;
    p(3) = c;
    
    s_params = p;
    options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
    estimate = fminsearch(@myfit2,s_params,options,tmp_f',tmp_amp);
    out(n).a = estimate(1);
    out(n).b = estimate(2);
    out(n).c = estimate(3);
    out(n).FWHM = c*(2*sqrt(2*log(2)));
    fitted_curve = a.*exp(-(((f*1e6)-b).^2)./(2*(c^2)));
    plot(f*1e6,fitted_curve,'r')
    hold on
    out(n).Q = out(n).b/out(n).FWHM;
    display(sprintf('The quality factor Q = %g',out(n).Q))
end
hold off
title('Single-Sided Amplitude Spectrum of y(t)','FontSize',16)
xlabel('Wavelength (um)','FontSize',14)
ylabel('|Y(f)|','FontSize',14)
axis([1.5 1.6 0 1.1*max(amp)])
filename = sprintf('%s/FT_Hz_Timesignal_fit.eps',foldername);
print(f2,'-depsc','-r200','-painters',filename);

foldername1 = 'Results';
filename = sprintf('%s/Pump_probe_parameters.txt',foldername1);
fileID = fopen(filename,'w');
name = {'Probe' 'Pump'};

for n = 1:2
    fprintf(fileID,sprintf('%s parameters:\n',name{n}));
    fprintf(fileID,sprintf('    Amplitude   = %g\n',out(n).a));
    fprintf(fileID,sprintf('    Wavelength  = %gum\n',out(n).b));
    fprintf(fileID,sprintf('    FWHM        = %gum\n',out(n).FWHM));
    fprintf(fileID,sprintf('    Quality factor Q    = %g\n',out(n).Q));
    if n == 1
        fprintf(fileID,sprintf('    Cavity photon lifetime  = %g s\n\n',out(n).Q*out(n).b*1e-6/(2*pi*si.c)));
    end
end
fprintf(fileID,sprintf('\nProbe/Pump amplitude = %g %%%%\n',100*out(1).a/out(2).a));
fprintf(fileID,sprintf('\nCriticial population level N2 = %g\n',av_E(3)));
fprintf(fileID,sprintf('Criticial inversion N2-N1 = %g\n',av_E(5)));
fclose(fileID);


end