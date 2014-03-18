% Gauss fit
clc
x = 1:1:5000;

actual.amp = 1;
actual.omega = 0.05;
actual.phase = 2.5;
actual.centre = 2000;
actual.width = 300;





actual.y = actual.amp.*(sin((actual.omega.*x)+actual.phase)).*exp(-((x-actual.centre).^2)./(2*(actual.width^2)));



% Now find estimates

% Amplitude
amp_min = min(actual.y);
amp_max = max(actual.y);

estimate.amp = (abs(amp_min)+abs(amp_max))/2;


% Frequency

sign_y = sign(actual.y);
zeros_y = sign_y(1:end-1)+sign_y(2:end);
ind = find(zeros_y == 0);
dif = ind(2:end)-ind(1:end-1);
lambda = 2*mean(dif);
estimate.omega = 2*pi/lambda;


% Phase
estimate.phase = 0;


% Centre

ind_min = find(actual.y == amp_min);
ind_max = find(actual.y == amp_max);
estimate.centre = (ind_min+ind_max)/2;


% Width
half_max = estimate.amp/2;

ind = find(actual.y >= half_max);
fwhm = ind(end)-ind(1);
estimate.width = fwhm/(2*sqrt(2*log(2)));

p(1) = estimate.amp;
p(2) = estimate.omega;
p(3) = estimate.phase;
p(4) = estimate.centre;
p(5) = estimate.width;
s_params = p;



options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
fitted = fminsearch(@myfit,s_params,options,x,actual.y);
format long

amp = fitted(1);
omega = fitted(2);
phase = fitted(3);
centre = fitted(4);
width = fitted(5);

fprintf('\nAmplitude = %g',amp)
fprintf('\nFrequency = %g',omega)
fprintf('\nPhase = %g',phase)
fprintf('\nCentre = %g',centre)
fprintf('\nWidth = %g\n\n',width)

estimate_y = amp.*(sin((omega.*x)+phase)).*exp(-((x-centre).^2)./(2*(width^2)));

figure(1)
plot(x,actual.y)
hold on
plot(x,estimate_y,'r')
hold off




