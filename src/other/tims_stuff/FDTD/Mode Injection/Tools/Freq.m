% Calculate frequencys
    c_si = 299792458;
    dt = 0.7;
    dx = 1.30537e-6/10/sqrt(11.68)
    dt_r = dt*dx/c_si;

    filename = 'Point_Hz.0.gpl';
    
    file = fopen(filename,'r');
    %[A,COUNTA] = textscan(file,'%n',1e5,'Delimiter','\n','CommentStyle','#');
    loop = 0;
    
    for a = 1:13
        tline = fgetl(file);
    end
    
    
    while 1
        loop = loop+1;
        tline = fgetl(file);
        if ~ischar(tline),   break, end
        
        if isempty(tline) == 0
            time_sig(loop,:) = str2num(tline);
        end
    end
    fclose(file);
    A   =   time_sig;
    x = A;
    N=length(x);
    Fs = 2;
    %this part of the code generates that frequency axis
    if mod(N,2)==0
        k=-N/2:N/2-1; % N even
    else
        k=-(N-1)/2:(N-1)/2; % N odd
    end
    T=N/Fs;
    freq=k/T;
    freq    =   freq/dt_r;
     y = x;
     zeroPadFactor = nextpow2(length(y)) + 3;
     [a,b] = positiveFFT_zero_padding(y,Fs,2^zeroPadFactor);
     figure(3)
     f = (Fs/dt_r/4)/2*linspace(0,1,NFFT*4);
     %b   =   (2*pi)./(dx./b);
     plot(c_si./f,abs(a))
    time_sig = time_sig(1:end);
    x= 1:length(time_sig);
    y = time_sig;
    x = x*dt_r*1e12;
    figure(1)
    plot(x,time_sig)
    yax = axis;
    axis([0 max(x) yax(3:end)])
    title('Hz field time signal','FontSize',16)
    xlabel('Time (ps)','FontSize',14)
    ylabel('Hz','FontSize',14)
    
    Fs = 1/dt_r;
    L = length(time_sig);
    NFFT = 2^nextpow2(L);
    Y = fft(y,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    f = c_si./f;
    amp = 2*abs(Y(1:NFFT/2+1));
    
    [pks,locs] = findpeaks(amp,'MINPEAKHEIGHT',max(amp)/10,'MINPEAKDISTANCE',100);
    figure(2)
    % Plot single-sided amplitude spectrum.
    plot(f*1e6,amp)
    hold on
    plot(f(locs)*1e6,pks*1.005,'k^','markerfacecolor',[1 0 0]);
    hold off
    for loop = 1:length(locs)
        lab = sprintf('%g',f(locs(loop))*1e6);
        text(f(locs(loop))*1e6+0.014,pks(loop),lab,'FontSize',14);
    end
    %lab = sprintf('%g',f(locs(2))*1e6);
    %text(f(locs(2))*1e6+0.014,pks(2),lab,'FontSize',14);    
    title('Single-Sided Amplitude Spectrum of y(t)','FontSize',16)
    xlabel('Wavelength (um)','FontSize',14)
    ylabel('|Y(f)|','FontSize',14)
    axis([1.2 1.7 0 1.1*max(2*abs(Y(1:NFFT/2+1)))])
    foldername = 'Results/Frequency Spectrum';
    if (exist(foldername, 'dir') == 0)
        mkdir(foldername)
    end
    
    filename = sprintf('%s/Hz_Timesignal.eps',foldername);
    print(f1,'-depsc','-r200','-painters',filename);
    filename = sprintf('%s/FT_Hz_Timesignal.eps',foldername);
    print(f2,'-depsc','-r200','-painters',filename);
    