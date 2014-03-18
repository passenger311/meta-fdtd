close all
clear all
clc

dx = 1e-8;
dt = dx*0.7/299792458;
t0 = 27838;

for loop = 1:10

tst = t0+(loop-1)*22;
tfi = tst+22;


fname{1} = sprintf('Hz_Slice_%g.0.gpl',tst);
fname{2} = sprintf('Hz_Slice_%g.0.gpl',tfi);

for i = 1:2
    tmp = f_readGPLFile2(fname{i});
    params = tmp{1};
    data = tmp{2};
    x = (1:length(data))*params(3)*dx;
    NFFT = 2^(nextpow2(length(x)) + 0);
    NFFT = length(x);
    tmp = fourier_trans(x,data');
    k = 2*pi*tmp(1,:);
    ftmp = fft(data,NFFT);
    ft = log(fft(data,NFFT));
    a{i} = real(ft);
    theta{i} = imag(ft);
end

Fs = 1/(x(2)-x(1));

f = Fs/2*linspace(0,2,length(ft));
t = dt*(tfi-tst);
rew = (theta{1}-theta{2})/t;
imw = (a{2}-a{1})/t;
st = floor(length(ft)/2);

fname2 = 'Clean_Complex_w__Mode_Tim_Leaky_lossy_500nm_3.txt';
tmp = dlmread(fname2,'\t');



figure(1)
plot(2*pi*f*1e-6,(rew)/2/pi/1e12,'ro')
hold on
plot(tmp(:,1),tmp(:,3))
xlim([0 2])
ylim([193.67 193.72])

figure(2)
plot(2*pi*f*1e-6,-2*imw/1e12,'bo')
hold on
plot(tmp(:,1),tmp(:,4))
xlim([0 2])
ylim([35 37])

figure(3)
plot(2*pi*f,2*abs(ftmp))
xlim([0 2e6])

end

figure
plot(x,data)