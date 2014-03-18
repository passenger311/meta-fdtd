close all
clear all

x       =   -8000:1:8000;

a1      =   1;
w1      =   0.101;
phi1    =   0;
b1      =   0;
c1      =   2500;

a2      =   1;
w2      =   0.1;
phi2    =   0;
b2      =   0;
c2      =   2500;

pulse1  =   a1.*sin((w1.*x)+phi1).*exp(-((x-b1).^2)./(2*(c1^2)));
pulse2  =   a2.*sin((w2.*x)+phi2).*exp(-((x-b2).^2)./(2*(c2^2)));
pulse3  =   pulse1+pulse2;
figure(1)
plot(x,pulse3)
Fs = 1;
y = pulse3;
N=length(x);
[YfreqDomain,freq] = centeredFFT(y,Fs);

figure(2)
plot(freq*2*pi,abs(YfreqDomain))
zeroPadFactor = nextpow2(length(y)) + 3;
[a,b] = positiveFFT_zero_padding(y,Fs,2^zeroPadFactor);
figure(3)
plot(b*2*pi,abs(a))
% p = unwrap(angle(X));
% p = fftshift(p);
% figure(3)
% plot(freq*2*pi,p)
return
p(1)    =   1;
p(2)    =   0.102;
p(3)    =   0;
p(4)    =   b1;
p(5)    =   c1;
% p(6)    =   a2;
% p(7)    =   w2;
% p(8)    =   0;
% p(9)    =   b2;
% p(10)   =   c2;



s_params = p;
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',5000,'Display','iter');
estimate = fminsearch(@myfit,s_params,options,x,pulse3);
n = 1;
%A_mat(n,:) = A;
a3          = estimate(1)
w3          = estimate(2)
phi3        = estimate(3)
b3          = estimate(4)
c3          = estimate(5)
% a4          = estimate(6)
% w4          = estimate(7)
% phi4        = estimate(8)
% b4          = estimate(9)
% c4          = estimate(10)

pulse4  =   a3.*sin((w3.*x)+phi3).*exp(-((x-b3).^2)./(2*(c3^2)));
%pulse4  =   pulse4+ a4.*sin((w4.*x)+phi4).*exp(-((x-b4).^2)./(2*(c4^2)));
figure(2)
plot(x,pulse3)
hold on
plot(x,pulse4,'r')
hold off

figure(3)
plot(x,(pulse4-pulse3))

y = pulse4-pulse3;

zeroPadFactor = nextpow2(length(y)) + 3;
[a,b] = positiveFFT_zero_padding(y,Fs,2^zeroPadFactor);
figure(4)
plot(b*2*pi,abs(a))
axis([0.05 0.15 0 max(abs(a))])