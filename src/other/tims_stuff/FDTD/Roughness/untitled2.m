close all
clear all
clc


dx = 0.01;
x = -20:dx:20;
L = length(x);
NFFT = 2^nextpow2(L);
k = 4;

y = cos(k*x);

figure
plot(x,y)

f = fft(y);

Fs = 1/dx;

N = length(f);
k = Fs/2*linspace(0,2,length(f));

figure
plot(2*pi*k,abs(f))

tmp = fourier_trans(x,y);

figure
plot(2*pi*tmp(1,:),abs(tmp(2,:)))