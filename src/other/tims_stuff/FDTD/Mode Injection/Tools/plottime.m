close all
clear all
clc

fname = 'Point_Hz.0.gpl';
fid = fopen(fname,'r');

for FLOOP = 1:13
    tline = fgetl(fid);
end
n = 0;
while ischar(tline)
    n = n + 1;
    tline = fgetl(fid);
    if ischar(tline)
        data(n) = str2num(tline);
    end
end
fclose(fid)


dx = 1e-8;
s = 0.7;
c = 299792458;
dt = 4*dx*s/c;


t = dt*(1:length(data))';
plot(t*1e12,data)

ind = find(abs(t-2e-12) == min(abs(t-2e-12)));

four = fourier_trans(t(ind:end)',data(ind:end));

figure(2)
semilogy(c./four(1,:),four(2,:))
xlim([1e-6 2e-6])
