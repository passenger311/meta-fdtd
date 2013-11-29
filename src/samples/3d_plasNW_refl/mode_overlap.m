function mode_overlap(file,file_mode,freqnumber,fig_on);
size_pml = 11;

fid1 = fopen('lambda2.in', 'r');
if (fid1==-1) error('file lambda.in does not exits'); end;
k=1;
while feof(fid1) ==0
   tline = fgetl(fid1);
   lambda(k)=sscanf(tline,'%e');
   k=k+1;
end
fclose(fid1);

fid1 = fopen('data.save', 'r');
if (fid1==-1) error('file data.save does not exits'); end;
tline = fgetl(fid1);
dx=sscanf(tline,'%e');
fclose(fid1);

%%% Read complex mode profiles for E and H from mode solver

% Read real part of E mode profile
filename=sprintf('Emode%s_re.set',file_mode);
%fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange1(:,:)=reshape(range,3,3)';
clear tmp;

tmp2 = textscan(fid, '%n %n %n ');
fclose(fid);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2);
   tmp(:,j) = tmp2{j};
end

x1 = (facerange1(1,1):facerange1(1,3):facerange1(1,2)).*dx;
y1 = (facerange1(2,1):facerange1(2,3):facerange1(2,2)).*dx;
z1 = (facerange1(3,1):facerange1(3,3):facerange1(3,2)).*dx;
if size(x1,2)==1
  dim11=y1;
  dim21=z1;
end
if size(y1,2)==1
  dim11=x1;
  dim21=z1;
end
if size(z1,2)==1
  dim11=x1;
  dim21=y1;
end
Emode = reshape(tmp,size(dim11,2),size(dim21,2),3);


% Read imaginary part of E mode profile
filename=sprintf('Emode%s_im.set',file_mode);
%fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange1(:,:)=reshape(range,3,3)';
clear tmp;

tmp2 = textscan(fid, '%n %n %n ');
fclose(fid);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2);
   tmp(:,j) = tmp2{j};
end

Emode = Emode+i*reshape(tmp,size(dim11,2),size(dim21,2),3);
Emode(:,:,3)=-Emode(:,:,3); %maybe the sign of Ez is calculated incorrectly in the mode solver? (same for Hz)


% Read real part of H mode profile
filename=sprintf('Hmode%s_re.set',file_mode);
%fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange1(:,:)=reshape(range,3,3)';
clear tmp;

tmp2 = textscan(fid, '%n %n %n ');
fclose(fid);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2);
   tmp(:,j) = tmp2{j};
end

Hmode = reshape(tmp,size(dim11,2),size(dim21,2),3);


% Read imaginary part of H mode profile
filename=sprintf('Hmode%s_im.set',file_mode);
%fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange1(:,:)=reshape(range,3,3)';
clear tmp;

tmp2 = textscan(fid, '%n %n %n ');
fclose(fid);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2);
   tmp(:,j) = tmp2{j};
end

Hmode = Hmode+i*reshape(tmp,size(dim11,2),size(dim21,2),3);
Hmode(:,:,3)=-Hmode(:,:,3);


%%% Read in diagmode result of FDTD
filename=sprintf('F_%s_%i.set',file,freqnumber);
%fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange(:,:)=reshape(range,3,3)';
clear tmp;

tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n %n ');
fclose(fid);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2);
   tmp(:,j) = tmp2{j};
end
E = zeros(size(tmp,1),3);
H = E;
E(:,1) = tmp(:,1)+i*tmp(:,2);
E(:,2) = tmp(:,3)+i*tmp(:,4);
E(:,3) = tmp(:,5)+i*tmp(:,6);
H(:,1) = tmp(:,7)+i*tmp(:,8);
H(:,2) = tmp(:,9)+i*tmp(:,10);
H(:,3) = tmp(:,11)+i*tmp(:,12);

x = (facerange(1,1):facerange(1,3):facerange(1,2)).*dx;
y = (facerange(2,1):facerange(2,3):facerange(2,2)).*dx;
z = (facerange(3,1):facerange(3,3):facerange(3,2)).*dx;
if size(x,2)==1
  dim1=y;
  dim2=z;
end
if size(y,2)==1
  dim1=x;
  dim2=z;
end
if size(z,2)==1
  dim1=x;
  dim2=y;
end

Ef = reshape(E,size(dim1,2),size(dim2,2),3);
Hf = reshape(H,size(dim1,2),size(dim2,2),3);

%{
% bringing zero phase to position of E-field maximum of Efmode
maxEf=max(max(Ef(1:floor(end/2)+1,:,1:2),[],2),[],1);
anEf=angle(maxEf);
[tmp,max_comp]=max(maxEf);
Efmode = Ef.*exp(-i*anEf(max_comp));
Hfmode = Hf.*exp(-i*anEf(max_comp));
%}

% defining fields for backward mode propagation
Hmode_b=-Hmode;
Emode_b=Emode;
Hmode_b(:,:,3)=-Hmode_b(:,:,3);
Emode_b(:,:,3)=-Emode_b(:,:,3);

% averaging mode profile from mode solver to Yee-cube center (as this is done in FDTD results)
Emode(:,:,1)=0.5*(Emode(:,:,1)+circshift(Emode(:,:,1),[1 0]));
Emode(:,:,2)=0.5*(Emode(:,:,2)+circshift(Emode(:,:,2),[0 1]));
Hmode(:,:,1)=0.5*(Hmode(:,:,1)+circshift(Hmode(:,:,1),[0 1]));
Hmode(:,:,2)=0.5*(Hmode(:,:,2)+circshift(Hmode(:,:,2),[1 0]));

Emode_b(:,:,1)=0.5*(Emode_b(:,:,1)+circshift(Emode_b(:,:,1),[1 0]));
Emode_b(:,:,2)=0.5*(Emode_b(:,:,2)+circshift(Emode_b(:,:,2),[0 1]));
Hmode_b(:,:,1)=0.5*(Hmode_b(:,:,1)+circshift(Hmode_b(:,:,1),[0 1]));
Hmode_b(:,:,2)=0.5*(Hmode_b(:,:,2)+circshift(Hmode_b(:,:,2),[1 0]));

% making Emode and Hmode complex (for function max to give desired results)
Emode = Emode + i*1e-100;
Hmode = Hmode + i*1e-100;

% making Emode_b and Hmode_b complex (for functions max to give desired results)
Emode_b = Emode_b + i*1e-100;
Hmode_b = Hmode_b + i*1e-100;

%{
% bringing zero phase back to position of field maximum of Emode
maxEf=max(max(Emode(:,:,1:2),[],2),[],1);
anEf=angle(maxEf);
[tmp,max_comp]=max(maxEf);
Emode = Emode.*exp(-i*anEf(max_comp));
Hmode = Hmode.*exp(-i*anEf(max_comp));

% making Emode and Hmode complex (for functions max to give desired results)
Emode = Emode - i*1e-100;
Hmode = Hmode - i*1e-100;

% bringing zero phase back to position of field maximum of Emode_b
maxEf=max(max(Emode_b(:,:,1:2),[],2),[],1);
anEf=angle(maxEf);
[tmp,max_comp]=max(maxEf);
Emode_b = Emode_b.*exp(-i*anEf(max_comp));
Hmode_b = Hmode_b.*exp(-i*anEf(max_comp));

% making Emode_b and Hmode_b complex (for functions max to give desired results)
Emode_b = Emode_b - i*1e-100;
Hmode_b = Hmode_b - i*1e-100;
%}


% Calculate cross poynting flux of FDTD and mode solver results for forward mode
overlapEmode=cross(Emode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:),conj(Hmode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:)),3);
overlapEfmode=cross(Efmode,conj(Hfmode),3);
norm_Emode=abs(sqrt(real(sum(sum(overlapEmode(:,:,3))))));
norm_Efmode=abs(sqrt(real(sum(sum(overlapEfmode(:,:,3))))));
fprintf('\nOverall power flow FDTD (+1 fw, -1 bw): %7.3d\n', real(sum(sum(overlapEfmode(:,:,3))))/2);
prompt = 'Use a larger value for the power flux?';
norm_Efmode=max(abs(sqrt(input(prompt))),norm_Efmode);

%{
% Loading or writing external normalisation constant
filename = sprintf('../%s/NormS_%s_%i.in',file_mode,file,freqnumber);
if exist(filename, 'file') == 2
  fprintf('Loading from file %s\n',filename)
  fid = fopen(filename, 'r');
  tline = fgetl(fid1);
  norm_Efmode=sscanf(tline,'%e');
  fprintf('Loaded power flow: %7.3d\n', norm_Efmode);
else
  fprintf('Writing file %s\n',filename)
  fid = fopen(filename, 'w');
  fprintf(fid,'%f\n',norm_Efmode);
end
fclose(fid1);
%}

overlap=cross(Emode(size_pml:end-size_pml+1,size_pml:size(Emode,2)-size_pml+1,:),conj(Hfmode),3);
ol1=sum(sum(overlap(:,:,3)));
overlap=cross(conj(Efmode),Hmode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:),3);
ol2=sum(sum(overlap(:,:,3)));
fprintf('\nForward mode overlap ave  = %7.3f %% at angle %6.1f degrees\n', abs(0.5*(ol1+ol2)/(norm_Emode*norm_Efmode))^2*100,angle(ol1+ol2)*180/pi);

% Calculate cross poynting flux of FDTD and mode solver results for backward mode
overlapEmode_b=cross(Emode_b(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:),conj(Hmode_b(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:)),3);
norm_Emode_b=abs(sqrt(real(sum(sum(overlapEmode_b(:,:,3))))));
overlap=cross(Emode_b(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:),conj(Hfmode),3);
ol1_b=sum(sum(overlap(:,:,3)));
overlap=cross(conj(Efmode),Hmode_b(size_pml:end-size_pml+1,size_pml:end-size_pml+1,:),3);
ol2_b=sum(sum(overlap(:,:,3)));
fprintf('Backward mode overlap ave = %7.3f %% at angle %6.1f degrees\n', abs(0.5*(ol1_b+ol2_b)/(norm_Emode_b*norm_Efmode))^2*100, angle(ol1_b+ol2_b)*180/pi);
fprintf('-------------------------------------------------------------\n\n')

if (fig_on==1); 
  % Plot comparison of real parts of Ex,Ey and imaginary part of Ez field profiles

  figure(1)
  subplot(3,2,1);
  pcolor(real(Emode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,1)'))
  title('Mode solver','FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,2);
  pcolor(real(Efmode(:,:,1)'))
  title(sprintf('FDTD (overlap %3.1f%%)',abs(0.5*(ol1+ol2)/(norm_Emode*norm_Efmode))^2*100),'FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,3);
  pcolor(real(Emode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,2)'))
  shading interp;
  subplot(3,2,4);
  pcolor(real(Efmode(:,:,2)'))
  shading interp;
  subplot(3,2,5);
  pcolor(imag(Emode(size_pml:end-size_pml+1,size_pml:end-size_pml+1,3)'))
  shading interp;
  subplot(3,2,6);
  pcolor(imag(Efmode(:,:,3)'))
  shading interp;
  filename = sprintf('E_%s_%s_%4.1fnm.png',file,file_mode,lambda(freqnumber));
  print(gcf, '-dpng', '-r100', filename);

  % Plot comparison of real parts of Hx,Hy and imaginary part of Hz field profiles  figure(2)
  subplot(3,2,1);
  pcolor(real(Hmode(size_pml:size(Emode,1)-size_pml+1,size_pml:size(Emode,2)-size_pml+1,1)'))
  title('Mode solver','FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,2);
  pcolor(real(Hfmode(:,:,1)'))
  title(sprintf('FDTD (overlap %3.1f%%)',abs(0.5*(ol1+ol2)/(norm_Emode*norm_Efmode))^2*100),'FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,3);
  pcolor(real(Hmode(size_pml:size(Emode,1)-size_pml+1,size_pml:size(Emode,2)-size_pml+1,2)'))
  shading interp;
  subplot(3,2,4);
  pcolor(real(Hfmode(:,:,2)'))
  shading interp;
  subplot(3,2,5);
  pcolor(imag(Hmode(size_pml:size(Emode,1)-size_pml+1,size_pml:size(Emode,2)-size_pml+1,3)'))
  shading interp;
  subplot(3,2,6);
  pcolor(imag(Hfmode(:,:,3)'))
  shading interp;

  filename = sprintf('H_%s_%s_%4.1fnm.png',file,file_mode,lambda(freqnumber));
  print(gcf, '-dpng', '-r100', filename);

  % Plot real part of Ex,Ey,Hx,Hy and imaginary part of Ez and Hz field profiles from FDTD
  figure(3)
  subplot(3,2,1);
  pcolor(real(Efmode(:,:,1)'))
  title('FDTD E','FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,2);
  pcolor(real(Hfmode(:,:,1)'))
  title('FDTD H','FontWeight','bold','FontSize',12)
  shading interp;
  subplot(3,2,3);
  pcolor(real(Efmode(:,:,2)'))
  shading interp;
  subplot(3,2,4);
  pcolor(real(Hfmode(:,:,2)'))
  shading interp;
  subplot(3,2,5);
  pcolor(imag(Efmode(:,:,3)'))
  shading interp;
  subplot(3,2,6);
  pcolor(imag(Hfmode(:,:,3)'))
  shading interp;

  filename = sprintf('Mode_%s_%4.1fnm.png',file,lambda(freqnumber));
  print(gcf, '-dpng', '-r100', filename);

end

