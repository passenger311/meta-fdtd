% Read file names
close all
clear all
clc
tstart = tic;
% Read params

filename    =   'sim_params.txt';
fileID      =   fopen(filename,'r');
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
typ         =   in{1}{1};
typ         =   isequal(typ,'dx');
if typ == 1
    resolution  =   str2double(in{1}{2});
else
    res     =   str2double(in{1}{2});
end

tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xstart      =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xfinish     =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
central_lambda  =   str2double(in{1}{2});
fclose(fileID);


lightSI = 299792458;
%central_omega       =   1.4e15;
%central_freq        =   central_omega/(2*pi);
central_lambda      =   central_lambda*1e-6;

if  typ == 1
    resolution  =   resolution*1e-6;
else
    
    resolution          = central_lambda/res/(sqrt(11.68));
end
Courant             = 0.7;
TimeStepREAL        = (Courant*resolution)/lightSI;
time                = 0;
time                = round(time/TimeStepREAL);
k0 = 2*pi/central_lambda;
% cd 28388.localhost/Data/
list = dir('Hz_*.gpl');
num = max(size(list));

for n = 1:num
    name = list(n).name;
    ind1 = strfind(name,'_');
    ind2 = strfind(name,'.0');
    str_num = str2num(name(ind1+1:ind2-1));
    num_list(n,1) = {str_num};

    num_list(n,2) = {(name(ind1+1:ind2-1))};
end

num_list = sortrows(num_list,1);
ind = [num_list{:,1}] > time;
num_list = num_list(ind,:);

num = max(size(num_list));
per = 0;
fprintf('Fitting ... %0.1f%%\n',per);

for n = 1:3
    clc
    pern     =   floor(1000*n/num);
    pern     =   pern/10;
    
    if pern>per
        
        fprintf('Fitting ... %0.1f%%\n',pern);
        per = pern;
        
    end
    
    ncyc = num_list{n,2};
    filename = sprintf('Hz_%s.0.gpl',ncyc);
    
    file = fopen(filename,'r');
    [A,COUNTA] = textscan(file,'%n',1e5,'Delimiter','\n','CommentStyle','#');
    fclose(file);
    A = A{1};
    
    ind  = isnan(A) == 0;
    A = A(ind);
    A_full = A;
    A = A(xstart:xfinish);
    A_mat{n} = A;
end
C = xstart:1:xfinish;
figure(1)
plot(C',A);
c1  =  A_mat{1};
c2  =   A_mat{2};
a1  =   fft(c1);
a2  =   fft(c2);
cr  =   a1.*conj(a2);
ouput   =   fft(cr);
a1 = zeros(1,1000);
a2  =   a1;
a1(101:200)     =   ones(1,100);
a2(401:500)     =   ones(1,100);
c1  =   fft(a1);
c2  =   fft(a2);
cr  =   c1.*conj(c2);
plot(abs(cr))


%autocorr([a1' a2']);
xcorr([a1 a2])