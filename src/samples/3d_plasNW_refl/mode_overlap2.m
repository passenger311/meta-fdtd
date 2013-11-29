function mode_overlap2(file_mode,file_mode2)

fid1 = fopen('data.save', 'r');
if (fid1==-1) error('file data.save does not exits'); end;
tline = fgetl(fid1);
dx=sscanf(tline,'%e');
fclose(fid1);

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

Emode2 = reshape(tmp,size(dim11,2),size(dim21,2),3);
Emode = Emode+i*Emode2;


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

Hmode2 = reshape(tmp,size(dim11,2),size(dim21,2),3);
Hmode= Hmode+i*Hmode2;
%Hmode = Hmode*exp(i*pi/2);
%Emode = Emode*exp(i*pi/2);


filename=sprintf('Emode%s_re.set',file_mode2);
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

Efmode = reshape(tmp,size(dim11,2),size(dim21,2),3);

filename=sprintf('Emode%s_im.set',file_mode2);
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

Efmode2 = reshape(tmp,size(dim11,2),size(dim21,2),3);
Efmode = Efmode+i*Efmode2;


filename=sprintf('Hmode%s_re.set',file_mode2);
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

Hfmode = reshape(tmp,size(dim11,2),size(dim21,2),3);

filename=sprintf('Hmode%s_im.set',file_mode2);
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

Hfmode2 = reshape(tmp,size(dim11,2),size(dim21,2),3);
Hfmode = Hfmode+i*Hfmode2;

%{
overlap=Emode.*conj(Efmode);
norm_Emode=sqrt(sum(sum(sum(abs(Emode).^2))));
norm_Efmode=sqrt(sum(sum(sum(abs(Efmode).^2))));
fprintf('\nNorm Emode = %4.2f\n', norm_Emode)
fprintf('Norm Efmode = %4.2f\n', norm_Efmode)
fprintf('  Mode overlap E = %4.1f per cent\n', abs(sum(sum(sum(overlap)))/(norm_Emode*norm_Efmode)*100))

overlap=Hmode.*conj(Hfmode);
norm_Hmode=sqrt(sum(sum(sum(abs(Hmode).^2))));
norm_Hfmode=sqrt(sum(sum(sum(abs(Hfmode).^2))));
fprintf('\nNorm Hmode = %4.2f\n', norm_Hmode)
fprintf('Norm Hfmode = %4.2f\n', norm_Hfmode)
fprintf('  Mode overlap H = %4.1f per cent\n', abs(sum(sum(sum(overlap)))/(norm_Hmode*norm_Hfmode)*100))
%}
overlap = cross(Emode,conj(Hfmode),3);
overlap2 = cross(conj(Efmode),Hmode,3);
overlapEmode=cross(Emode,conj(Hmode),3);
overlapEfmode=cross(Efmode,conj(Hfmode),3);
ol1=sum(sum(overlap(:,:,3)));
ol2=sum(sum(overlap2(:,:,3)));
norm_Emode=abs(sqrt(real(sum(sum(overlapEmode(:,:,3))))));
norm_Efmode=abs(sqrt(real(sum(sum(overlapEfmode(:,:,3))))));

%fprintf('\nNorm Emode = %4.2f\n', norm_Emode)
%fprintf('Norm Efmode = %4.2f\n', norm_Efmode)
fprintf('\nMode overlap = %7.3f %% at angle %6.1f degrees\n', abs(0.5*(ol1+ol2)/(norm_Emode*norm_Efmode))^2*100,angle(ol1+ol2)*180/pi);
%fprintf('\nMode overlap = %4.1f per cent\n', sum(sum(overlap(:,:,3)))/(norm_Emode*norm_Efmode)*100)
%overlap=cross(Efmode,conj(Hmode),3);
%fprintf('Mode overlap = %4.1f per cent\n', sum(sum(overlap(:,:,3)))/(norm_Emode*norm_Efmode)*100)

