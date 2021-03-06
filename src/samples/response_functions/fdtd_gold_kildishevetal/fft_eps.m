clear all;
fid = fopen('data.save','r');
tline = fgetl(fid); d_mat = sscanf(tline,'%e');
tline = fgetl(fid); n_bg = sscanf(tline,'%e');
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

filename=sprintf('fft_r.pspec');
fprintf(1,'\nloading file: %s', filename);
fid1 = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

for j=1:7
  tline = fgetl(fid1);
end

tmp2 = textscan(fid1, '%n %n %n %n %n');
fclose(fid1);

fft_r = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   fft_r(:,j) = tmp2{j};
end
r = zeros(size(fft_r,1),1);
rp = zeros(size(fft_r,1),1);
r(:) = fft_r(:,2).*exp(i.*fft_r(:,3));


filename=sprintf('fft_t.pspec');
fprintf(1,'\nloading file: %s', filename);
fid1 = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

for j=1:7
  tline = fgetl(fid1);
end

tmp2 = textscan(fid1, '%n %n %n %n %n');
fclose(fid1);

fft_t = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   fft_t(:,j) = tmp2{j};
end
t = zeros(size(fft_t,1),1);
tp= zeros(size(fft_t,1),1);
t(:) = fft_t(:,2).*exp(i*fft_t(:,3));
fft1(:,1) = fft_t(:,1);

rp(:) = r(:).*exp(-i.*2.*pi.*n_bg.*fft1(:,1).*4);
tp(:) = t(:).*exp(-i.*2.*pi.*n_bg.*fft1(:,1).*2);

eps_th = eps_inf - o_D^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*g_D.*2.*pi.*fft1(:,1)) - de_L1.*o_L1^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*2.*g_L1.*2.*pi.*fft1(:,1)-o_L1^2);% -  de_L2.*o_L2^2./(2.*pi.*fft1(:,1).*2.*pi.*fft1(:,1)-i.*2.*g_L2.*2.*pi.*fft1(:,1)-o_L2^2);

eps_th = conj(eps_th);
n_th2 = sqrt(eps_th);

n = ( acos ( ( 1- ( rp(:).*rp(:) - tp(:).*tp(:) ) ) ./ ( 2.*tp(:)) ) + 2.*pi.*0 ) ./ ( 2.*pi.*fft1(:,1).*d_mat);

eps=n.*n;

figure(1)
plot(dx./fft1(:,1),real(eps),'o',dx./fft1(:,1),real(eps_th));
axis([350 800 -10 10])
axis 'auto y'
figure(2)
plot(dx./fft1(:,1),imag(eps),'o',dx./fft1(:,1),imag(eps_th));
axis([350 800 -10 10])
axis 'auto y'
figure(3)
plot(dx./fft1(:,1),real(n),'o',dx./fft1(:,1),real(n_th2));
axis([350 800 -.5 0]);
axis 'auto y'
figure(4)
plot(dx./fft1(:,1),imag(n),'o',dx./fft1(:,1),imag(n_th2));
axis([350 800 -10 10])
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

