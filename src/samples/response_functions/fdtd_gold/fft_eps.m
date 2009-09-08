clear all;
fid = fopen('data.save','r');
tline = fgetl(fid); dx = sscanf(tline,'%e');
tline = fgetl(fid); invl = sscanf(tline,'%e');
tline = fgetl(fid); eps_inf = sscanf(tline,'%e');
tline = fgetl(fid); o_D = sscanf(tline,'%e'); o_D = 2 * pi * o_D;
tline = fgetl(fid); g_D = sscanf(tline,'%e');
tline = fgetl(fid); de_L1 = sscanf(tline,'%e');
tline = fgetl(fid); o_L1 = sscanf(tline,'%e'); o_L1 = 2 * pi * o_L1;
tline = fgetl(fid); g_L1 = sscanf(tline,'%e');
%tline = fgetl(fid); de_L2 = sscanf(tline,'%e');
%tline = fgetl(fid); o_L2 = sscanf(tline,'%e'); o_L2 = 2 * pi * o_L2;
%tline = fgetl(fid); g_L2 = sscanf(tline,'%e');
fclose(fid);

filename=sprintf('fft_2.pspec');
fprintf(1,'\nloading file: %s', filename);
fid1 = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

for j=1:7
  tline = fgetl(fid1);
end

tmp2 = textscan(fid1, '%n %n %n %n %n');
fclose(fid1);

fft1 = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   fft1(:,j) = tmp2{j};
end
eps_th = eps_inf - o_D^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*g_D.*2.*pi.*fft1(:,1)) - de_L1.*o_L1^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*2.*g_L1.*2.*pi.*fft1(:,1)-o_L1^2);% -  de_L2.*o_L2^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*2.*g_L2.*2.*pi.*fft1(:,1)-o_L2^2);

eps_th = conj(eps_th);
n_th2 = sqrt(eps_th);

n=(log(fft1(:,2))./i+fft1(:,3))./(2.*pi.*fft1(:,1).*4);
eps=n.*n;

figure(1)
plot(dx./fft1(:,1),real(eps),dx./fft1(:,1),real(eps_th),'o');
axis([300 800 -10 10])
axis 'auto y'
figure(2)
plot(dx./fft1(:,1),imag(eps),dx./fft1(:,1),imag(eps_th),'o');
axis([300 800 -10 10])
axis 'auto y'
figure(3)
plot(dx./fft1(:,1),real(n),dx./fft1(:,1),real(n_th2),'o');
axis([300 800 -.5 0]);
axis 'auto y'
figure(4)
plot(dx./fft1(:,1),imag(n),dx./fft1(:,1),imag(n_th2),'o');
axis([300 800 -10 10])
axis 'auto y'

filename=sprintf('eps_real.dat');
fprintf(1,'\nwriting file: %s', filename);
fid1 = fopen(filename, 'w');
filename=sprintf('eps_imag.dat');
fprintf(1,'\nwriting file: %s\n', filename);
fid2 = fopen(filename, 'w');

for i=1:size(n,1)
  if ((dx/fft1(i,1)>=250) & (dx/fft1(i,1)<=1000))
    fprintf(fid1,'%e %e\n',dx/fft1(i,1),real(eps(i)));
    fprintf(fid2,'%e %e\n',dx/fft1(i,1),imag(eps(i)));
  end
end

fclose(fid1);
fclose(fid2);

