% This program aims to extract the complex frequency values from                          
% Field Vs Space plots. 

global x
global y
% Constants
lightSI = 299792458;                                              % [m/s] Speed of light in vacuum 

% Simulation parameters - Computational
Geometrylength        = 1024;                                    % [dx] The length of the geometry i.e. number of points which will appear in read in files
xstart                = 12;                                        % Lowest x-axis value of snapshot range   
xfinish               = 1710;                                     % Highest x-axis value of snapshot range
res                    = 10 ;                                   %[m]
resolution            = 1.55e-6/res/(sqrt(11.68))
frequencycentral      = lightSI/1.55e-6;
WavelengthcentralfreqREAL = (lightSI/frequencycentral);
Wavelengthcentralfreq = WavelengthcentralfreqREAL/resolution; % 113.89897;    % [dx] Gives the wavelength [grid units] of central frequency                                   
                               % [m] Wavelength of central frequency
Courant               = 0.7;                                      % Courant value used in simulation
TimeSeparation        = 40;                                    % [dt] Separation in time between snap shots 
maxfrac               = 0.1;                                     % Fraction of maxamplitude to use for maximum and minimum of range  
fprintf('Wavelengthcentralfreq %12.9E\n',Wavelengthcentralfreq);
%Converting units to strings so they can be put in file title
stxstart = sprintf('%d',xstart);
stxfinish = sprintf('%d',xfinish);

% Conversion - Converting from computational to SI units


% Conversion - Calculation
SpaceStepREAL = WavelengthcentralfreqREAL/Wavelengthcentralfreq;  % [m] Spatial Length of cell in SI units
TimeStepREAL = (Courant*SpaceStepREAL)/lightSI;                   % [s] Length of each time step in SI units
TimeSeparationREAL = TimeStepREAL*TimeSeparation;                 % [s] Separation in time between snap shots

% Constants within program
MaxWavevector = (2*pi)/SpaceStepREAL;                             % [rad/m] Maximum value of wavevector 



% Controlling input
st1='~/PhD/FDTD_projects/Test/Hz_6610.0.gpl';           % Snapshot at earliest time
st2='~/PhD/FDTD_projects/Test/Hz_6610.0.gpl';           % Snapshot at later time
                                         
numb1 = regexp(st1,'/','split');                                         % Split up pathway st1 using "/" as divider
numb2 = regexp(char(numb1(length(numb1))),'_','split');                  % Split up final element of st1 using "_" as divider
numb3 = regexp(char(numb2(length(numb2))),'\.','split');                 % Split up final element of numb2 using "." as divider 
stringst1 = char(numb3(1));                                              % Stringst1 = time step of results contained in file

numbA1 = regexp(st2,'/','split');
numbA2 = regexp(char(numbA1(length(numbA1))),'_','split');
numbA3 = regexp(char(numbA2(length(numbA2))),'\.','split');
stringst2 = char(numbA3(1));

xrange = xfinish-xstart+1;                                               % No. points in snapshot range
%nextpower2 = nextpow2(xrange)+3;                                        % Raises length of geometry to next power of two 
nextpower2 = 17;
N = 2^nextpower2;


% Controlling output
st3='~/PhD/FDTD_projects/Test/fourierSEPcheck-SPACE';
suffix = '-maximum';
ftype= '.gpl';
st4=strcat(st3,'-',stringst1,'-',stringst2,'T-',stxstart,'-',stxfinish,'X',suffix,ftype);
st3=strcat(st3,'-',stringst1,'-',stringst2,'T-',stxstart,'-',stxfinish,'X',ftype);

% Reading in output from simulations
file_1 = fopen(st1,'r');
[A,COUNTA] = textscan(file_1,'%n',N,'Delimiter','\n','CommentStyle','#');
fclose(file_1);

%File 2
file_2 = fopen(st2,'r');
[B,COUNTB] = textscan(file_2,'%n',N,'Delimiter','\n','CommentStyle','#');
fclose(file_2)


A=A{1};                     
B=B{1};
size(A)

A=A(xstart+2:xfinish-1); % Go from start+1 to remove the NAN which is a consequence of the empty line between header and data               
B=B(xstart+2:xfinish-1);
C = 1:1:length(A);
% figure(1);
% plot(C,A,C,B,'--')
% hold on
% a = (abs(min(A))+abs(max(A)))/2;
% A_s = sign(A);
% A_s2 = A_s(1:end-1)+A_s(2:end);
% ind = find(A_s2 == 0);
% ind2 = ind(2:end)-ind(1:end-1);
% ind2 = ind2(4:end-2);
% l_2 = mean(ind2);
% lambda = l_2*2;
% freq = 1/lambda;
% 
% w = 2*pi*freq;
% phi = 0;
% b = (find(A == max(A)) + find(A == min(A)))/2;
% ind = find(A >= a/2);
% fwhm = ind(end)-ind(1);
% c = fwhm/(2*sqrt(2*log(2)));
% p(1) = a;
% p(2) = w;
% p(3) = phi;
% p(4) = b;
% p(5) = c;
% s_params = p;
% options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',2000);
% estimate = fminsearch(@myfit,s_params,options,C',A);
% 
% x = C;
% y = A;
% a = estimate(1)
% w = estimate(2)
% phi = estimate(3)
% b = estimate(4)
% c = estimate(5)
% D = a.*sin((w.*x)+phi).*exp(-((x-b).^2)./(2*(c^2)));
% 
% plot(C,D,'r')
% hold off
% fileID = fopen('gaus_out.txt','w');
% for loop = 1:length(C)
%     fprintf(fileID,'%g\t%g\n',C(loop),A(loop));
% end
% fclose(fileID);
Fa=fft(A,2^nextpower2);           % Fourier transform the signal 
Fb=fft(B,2^nextpower2);    

Fa1 = Fa;                                           
Fb1 = Fb;

% Obtaining array of Wavevectors - Going from 0 - Max
Wavevector = ((0:((2^nextpower2)/2))/(2^nextpower2)).*MaxWavevector; % [rad/m] Range of Wavevectors extending from 0 to MaxWavevector


Fpa = 2*Fa1(1:(2^nextpower2)/2+1)/(2^nextpower2);                           % Multiply by 2 as you are only taking half of the points
Fpb = 2*Fb1(1:(2^nextpower2)/2+1)/(2^nextpower2);                           % Divide each by nextpower so you have magnitude per wavevector.

%figure(2);plot(Wavevector,abs(Fpa),Wavevector,abs(Fpb),'--')

absFpa=abs(Fpa);                                                            % Find the absolute value of all frequency component
maxamp=max(absFpa);                                                         % Find value of the largest absolute frequency component
range = maxamp*maxfrac;                                                     % Find minimum value of absolute frequency component to include
maxindex = find(absFpa==maxamp);                                            % Find the index the minimum absolute frequency component

% Array definitions which are required
phaseA = zeros(1,numel(Fpa));
phaseB = zeros(1,numel(Fpa));
Changephase = zeros(1,numel(Fpa));
ampA = zeros(1,numel(Fpa));
ampB = zeros(1,numel(Fpa));
OmegaR = zeros(1,numel(Fpa));
OmegaI = zeros(1,numel(Fpa));
Changeamp = zeros(1,numel(Fpa));

Frequency = zeros(1,numel(Fpa));
Nreal = zeros(1,numel(Fpa));
Nimag = zeros(1,numel(Fpa));

for u =2:numel(Fpa);
    phaseA(u) = atan2(imag(Fpa(u)),real(Fpa(u)));
    phaseB(u) = atan2(imag(Fpb(u)),real(Fpb(u)));
    Changephase(u) = phaseB(u)-phaseA(u);
    
    if Changephase(u)>pi;
                Changephase(u)=Changephase(u)-2*pi;
    elseif Changephase(u)<-pi;
                Changephase(u)=Changephase(u)+2*pi;
    end;
              
    ampA(u)=abs(Fpa(u));
    ampB(u)=abs(Fpb(u));
    Changeamp(u) = ampB(u)-ampA(u);
    
    OmegaR(u) = -1*Changephase(u)/(TimeSeparationREAL); % So A-B
    OmegaI(u) = (log(ampB(u)/ampA(u)))/(TimeSeparationREAL);
    
    
    % Calculating the real part of the refractive index
    
    Nreal(u) = Wavevector(u)/(Changephase(u)/(TimeSeparationREAL*lightSI)); % [Unitless] k/k0 
                                                                               % Calculate k0 by first calculating frequency and then converting
    %Nimag(u) = Changeamp(u)/((-1*Changephase(u)/(TimeSeparationREAL*lightSI))*lightSI*TimeSeparationREAL);
    
    Nimag(u) = (ampB(u)-ampA(u))/(Changephase(u));
    % Calculating frequency (ASSUMING phase velocity is equal to c
    Frequency(u) = (Wavevector(u)/(2*pi))*lightSI;
    
end;

out = fopen(st3,'wt');
fprintf(out,'# fourierspace.m\n')
fprintf(out,'# file_1 %s\n',st1');
fprintf(out,'# file_2 %s\n',st2');
fprintf(out,'# Computational Parameters\n');
fprintf(out,'# Geometrylength %12.6E Wavelengthcentralfreq %12.6E Courant %12.6E TimeSeparation %12.6E\n', Geometrylength,Wavelengthcentralfreq,Courant,TimeSeparation');
fprintf(out,'# SI Parameters\n');
fprintf(out,'# WavelengthcentralfreqREAL %12.6E SpaceStepREAL %12.6E TimeStepREAL %12.6E TimeSeparationREAL %12.6E xstart %12.6E xfinish %12.6E\n',WavelengthcentralfreqREAL, SpaceStepREAL, TimeStepREAL, TimeSeparationREAL, xstart, xfinish')

for u = 3:numel(Fpa);
    if abs(Fpa(u))>range;
%     fprintf(out,'Wavevector %12.6E ampA %12.6E ampB %12.6E ampA/ampB %18.12E log(ampA/ampB) %18.12E log(ampA(u)/ampB(u))/TimeSeparationREAL) %18.12E \n',Wavevector(u),ampA(u),ampB(u), ampA(u)/ampB(u), log(ampA(u)/ampB(u)),log(ampA(u)/ampB(u))/TimeSeparationREAL);  
%     fprintf(out,'Frequency %12.6E Changeamp %12.6E ampA %12.6E ampB %12.6E u %12.6E\n',Frequency(u),Changeamp(u),ampA(u),ampB(u),u); 
%     fprintf(out,'imag(Fpa(u)) %12.6E real(Fpa(u)) %12.6E abs(Fpa(u)) %12.6E\n',imag(Fpa(u)), real(Fpa(u)), abs(Fpa(u)))
%     fprintf(out,'real(Fpa) %12.6E imag(Fpb) %12.6E real(Fpa) %12.6E imag(Fpb) %12.6E\n',real(Fpa(u)), imag(Fpa(u)),real(Fpb(u)),imag(Fpb(u))); 
%     fprintf(out,'phaseA %12.6E phaseB %12.6E Changephase %12.6E ampA %12.6E ampB %12.6E\n',phaseA(u),phaseB(u),Changephase(u),ampA(u),ampB(u)); 
     fprintf(out,'Frequency(THz) %+12.6E Wavevector %+12.6E OmegaR %+12.6E OmegaI %+12.6E Nref %+12.6E Nreal %+12.6E Nimag %+12.6E ampA %+12.6E ampB %+12.6E Changephse %+12.6E\n',Frequency(u)/1.0E12,Wavevector(u),OmegaR(u),OmegaI(u),(Wavevector(u))/OmegaR(u),Nreal(u),Nimag(u),ampA(u), ampB(u), Changephase(u));   
    end;
end;
fclose(out)

out2 = fopen(st4,'wt');
fprintf(out2,'# fourierspace.m\n')
fprintf(out2,'# file_1 %s\n',st1');
fprintf(out2,'# file_2 %s\n',st2');
fprintf(out2,'# Computational Parameters\n');
fprintf(out2,'# Geometrylength %12.6E Wavelengthcentralfreq %12.6E Courant %12.6E TimeSeparation %12.6E\n', Geometrylength,Wavelengthcentralfreq,Courant,TimeSeparation');
fprintf(out2,'# SI Parameters\n');
fprintf(out2,'# WavelengthcentralfreqREAL %12.6E SpaceStepREAL %12.6E TimeStepREAL %12.6E TimeSeparationREAL %12.6E xstart %12.6E xfinish %12.6E\n',WavelengthcentralfreqREAL, SpaceStepREAL, TimeStepREAL, TimeSeparationREAL, xstart, xfinish')
fprintf(out2,'Frequency(THz) %12.6E Wavevector %12.6E OmegaR %12.6E OmegaI %12.6E Nref %12.6E Nreal %12.6E Nimag %12.6E ampA %12.6E ampB %12.6E\n',Frequency(maxindex)/1.0E12,Wavevector(maxindex),OmegaR(maxindex),OmegaI(maxindex),(Wavevector(maxindex))/OmegaR(maxindex),Nreal(maxindex),Nimag(maxindex),ampA(maxindex), ampB(maxindex));   
  
figure(2)
plot(Wavevector,ampA)
fclose(out2)
k0 = 2*pi/(WavelengthcentralfreqREAL)
ind = find(ampA == max(ampA))
k = Wavevector(ind)
Neff = k/k0


