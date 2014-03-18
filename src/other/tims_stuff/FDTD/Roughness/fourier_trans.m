function spec = fourier_trans(x,y)

zeroPadFactor = nextpow2(length(y)) + 0;

Fs      =   1/(x(2)-x(1));
L       =   length(y);
N       =   2^zeroPadFactor;
k       =   0:N-1;                          %create a vector from 0 to N-1
T       =   N/Fs;                           %get the frequency interval
X       =   fft(y,N)/L;                     % normalize the data

%only want the first half of the FFT, since it is redundant
cutOff  =   ceil(N/2);

%take only the first half of the spectrum
X       =   X(1:cutOff);

f       =   Fs/2*linspace(0,1,N/2+1);
f       =   f(1:end-1);

amp     =   (X(1:N/2));

spec = [f;amp];