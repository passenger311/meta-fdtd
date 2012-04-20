clear all;


fid = fopen('data.save','r');
tline = fgetl(fid); dx = sscanf(tline,'%e');
tline = fgetl(fid); d_mat = sscanf(tline,'%e');
tline = fgetl(fid); n_front = sscanf(tline,'%e');
tline = fgetl(fid); n_back = sscanf(tline,'%e');
tline = fgetl(fid); invl = sscanf(tline,'%e');
tline = fgetl(fid); d_tfsf_neg = sscanf(tline,'%e');
tline = fgetl(fid); d_tfsf_pos = sscanf(tline,'%e');
fclose(fid);

%% Load reference file

filename=sprintf('fft_ref.pspec');
fprintf(1,'\nloading file: %s', filename);
fid1 = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

for j=1:7
  tline = fgetl(fid1);
end

tmp2 = textscan(fid1, '%n %n %n %n %n');
fclose(fid1);

fft_ref = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   fft_ref(:,j) = tmp2{j};
end


%% Load non-normalized reflection file

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
r(:) = fft_r(:,2)./fft_ref(:,2).*exp(i.*(fft_r(:,3)-fft_ref(:,3)));


%% Load non-nomralized transmission file

filename=sprintf('fft_t.pspec');
fprintf(1,'\nloading file: %s\n', filename);
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
t(:) = fft_t(:,2)./fft_ref(:,2).*exp(i*(fft_t(:,3)-fft_ref(:,3)));
fft1(:,1) = fft_t(:,1);
lambda = dx./fft1(:,1);

%% Apply phase shift to complex reflection and transmission

rp(:) = r(:).*exp(-i.*2.*pi.*n_front.*fft1(:,1).*2.*d_tfsf_neg);
tp(:) = t(:).*exp(-i.*2.*pi.*fft1(:,1).*(n_back.*((d_tfsf_pos-2)-d_mat)+n_front.*(d_tfsf_neg-2)));

%% Parameter retrieval of refractive index

n = (acos((n_front-(n_front.*rp(:).*rp(:)-n_back.*tp(:).*tp(:))) ./((n_front+n_back).*tp(:)+(n_back-n_front).*tp(:).*rp(:)))-0*2*pi)./(2.*pi.*fft1(:,1).*d_mat);

%% Correct for phase jumps

%Choose wavelength region
[minimum, position1]=min(abs(dx./fft1(:,1)-800));
[minimum, position2]=min(abs(dx./fft1(:,1)-590));

%Force positive imaginary refractive index (not applicable for gain media)
n = n .* ((imag(n)>=0) - (imag(n)<0));

%Automatially adjust m factor 
for j=position1+1:size(n)
  while(true)
    if (real(n(j)-n(j-1))) > 1*pi./(2.*pi.*fft1(j,1).*d_mat)
      n(j) = n(j) - 2*pi./(2.*pi.*fft1(j,1).*d_mat);
    elseif (real(n(j)-n(j-1))) < -1*pi./(2.*pi.*fft1(j,1).*d_mat)
      n(j) = n(j) + 2*pi./(2.*pi.*fft1(j,1).*d_mat);
    else
      break;
    end
  end
end

%% Parameter retrieval of impedance

z2 = (((1+rp(:)).*(1+rp(:))-tp(:).*tp(:))./(n_front^2.*(1-rp(:)).*(1-rp(:))-n_back^2.*tp(:).*tp(:)));
z = sqrt(z2);

%Force positive real part of impedance
z = z .* ((real(z)>=0) - (real(z)<0));
%not necessarily applicable to gain media

n2 = 1./(2.*pi.*i.*fft1(:,1).*d_mat).*log(tp(:)./(1-rp(:).*(z(:)-1)./(z(:)+1)));
for j=position1+1:size(n)
  while(true)
    if (real(n2(j)-n2(j-1))) > 1*pi./(2.*pi.*fft1(j,1).*d_mat)
      n2(j) = n2(j) - 2*pi./(2.*pi.*fft1(j,1).*d_mat);
    elseif (real(n2(j)-n2(j-1))) < -1*pi./(2.*pi.*fft1(j,1).*d_mat)
      n2(j) = n2(j) + 2*pi./(2.*pi.*fft1(j,1).*d_mat);
    else
      break;
    end
  end
end

%Determine permittivity and permeablity
eps = n2./z;
mu = n2.*z;

%% Write data to files

%Open files
filename=sprintf('eps.dat');
fprintf(1,'\nwriting file: %s', filename);
fid1 = fopen(filename, 'w');
filename=sprintf('n.dat');
fprintf(1,'\nwriting file: %s', filename);
fid2 = fopen(filename, 'w');
filename=sprintf('n2.dat');
fprintf(1,'\nwriting file: %s', filename);
fid3 = fopen(filename, 'w');
filename=sprintf('FOM.dat');
fprintf(1,'\nwriting file: %s\n', filename);
fid4 = fopen(filename, 'w');
filename=sprintf('mu.dat');
fprintf(1,'\nwriting file: %s', filename);
fid5 = fopen(filename, 'w');
filename=sprintf('z.dat');
fprintf(1,'\nwriting file: %s', filename);
fid6 = fopen(filename, 'w');

filename=sprintf('t_abs.dat');
fprintf(1,'\nwriting file: %s', filename);
fid7 = fopen(filename, 'w');
filename=sprintf('t_angle.dat');
fprintf(1,'\nwriting file: %s\n', filename);
fid8 = fopen(filename, 'w');
filename=sprintf('r_abs.dat');
fprintf(1,'\nwriting file: %s', filename);
fid9 = fopen(filename, 'w');
filename=sprintf('r_angle.dat');
fprintf(1,'\nwriting file: %s\n', filename);
fid10 = fopen(filename, 'w');
filename=sprintf('abs.dat');
fprintf(1,'\nwriting file: %s\n', filename);
fid11 = fopen(filename, 'w');


fprintf(fid1,'%e %e\n',[lambda(position1:position2),real(eps(position1:position2))]');
fprintf(fid1,'\n');
fprintf(fid1,'%e %e\n',[lambda(position1:position2),imag(eps(position1:position2))]');

fprintf(fid2,'%e %e\n',[lambda(position1:position2),real(n(position1:position2))]');
fprintf(fid2,'\n');
fprintf(fid2,'%e %e\n',[lambda(position1:position2),imag(n(position1:position2))]');

fprintf(fid3,'%e %e\n',[lambda(position1:position2),real(n2(position1:position2))]');
fprintf(fid3,'\n');
fprintf(fid3,'%e %e\n',[lambda(position1:position2),imag(n2(position1:position2))]');

fprintf(fid4,'%e %e\n',[lambda(position1:position2),abs(real(n2(position1:position2))./imag(n2(position1:position2)))]');

fprintf(fid5,'%e %e\n',[lambda(position1:position2),real(mu(position1:position2))]');
fprintf(fid5,'\n');
fprintf(fid5,'%e %e\n',[lambda(position1:position2),imag(mu(position1:position2))]');

fprintf(fid6,'%e %e\n',[lambda(position1:position2),real(z(position1:position2))]');
fprintf(fid6,'\n');
fprintf(fid6,'%e %e\n',[lambda(position1:position2),imag(z(position1:position2))]');

fprintf(fid7,'%e %e\n',[lambda(position1:position2),abs(tp(position1:position2)).*abs(tp(position1:position2))]');
fprintf(fid8,'%e %e\n',[lambda(position1:position2),angle(tp(position1:position2))]');
fprintf(fid9,'%e %e\n',[lambda(position1:position2),abs(rp(position1:position2)).*abs(rp(position1:position2))]');
fprintf(fid10,'%e %e\n',[lambda(position1:position2),angle(rp(position1:position2))]');

absorption = 1.0-abs(rp).*abs(rp)-abs(tp).*abs(tp);
fprintf(fid11,'%e %e\n',[lambda(position1:position2),absorption(position1:position2)]');

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid10);
fclose(fid11);

exit
