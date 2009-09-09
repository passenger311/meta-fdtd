function vismode(freqnumber,fignum);

fid = fopen('data.save','r');
tline = fgetl(fid);
dx = sscanf(tline, '%e');
fclose(fid);

fid1 = fopen('lambda2.in', 'r');
if (fid1==-1) error('file lambda2.in does not exits'); end;
k=1;
while feof(fid1) ==0
   tline = fgetl(fid1);
   lambda(k)=sscanf(tline,'%e');
   k=k+1;
end
fclose(fid1);


filename=sprintf('EHT_%i.set',freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;
filename=sprintf('../EHT_%i.set',freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid1 = fopen(filename, 'r');
if (fid1==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
tmp = textscan(fid1, '%s',1); clear tmp;

tmp = textscan(fid1, '%n %n %n %n %n %n %n %n %n', 1);
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2)
   range(j) = tmp{j};
end
facerange(:,:)=reshape(range,3,3)';
clear tmp; clear tmp1;
tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n');
tmp3 = textscan(fid1, '%n %n %n %n %n %n %n %n');
fclose(fid); fclose(fid1);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
tmp1 = zeros(size(tmp3{1},1),size(tmp3,2));
for j = 1 : size(tmp2,2)
   tmp(:,j) = tmp2{j};
   tmp1(:,j) = tmp3{j};
end

Ex = tmp(:,1)+i*tmp(:,2);
Ey = tmp(:,3)+i*tmp(:,4);
Hx = tmp(:,5)+i*tmp(:,6);
Hy = tmp(:,7)+i*tmp(:,8);
Erefx = tmp1(:,1)+i*tmp1(:,2);
Erefy = tmp1(:,3)+i*tmp1(:,4);
Hrefx = tmp1(:,5)+i*tmp1(:,6);
Hrefy = tmp1(:,7)+i*tmp1(:,8);

Erx = sum(Erefx)/size(Erefx,1);
Ery = sum(Erefy)/size(Erefy,1);
Efx = reshape(Ex(:,1),floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));
Efy = reshape(Ey(:,1),floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));


x = dx*(facerange(1,1):facerange(1,3):facerange(1,2));
y = dx*(facerange(2,1):facerange(2,3):facerange(2,2));
figure(fignum);
surf(y,x,sqrt(abs(Efx.*Efx)+abs(Efy.*Efy))/sqrt(abs(Erx)^2+abs(Ery)^2));
a=get(gca,'zlim');
set(gca,'clim',[0,ceil(a(2))]);
title(sprintf('Mode wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',20);
shading interp;
axis tight;
colorbar;
view(90,-90);
xlabel('y in nm','FontWeight','bold','FontSize',16);
ylabel('x in nm','FontWeight','bold','FontSize',16);
set(gca,'DataAspectRatio',[1 1 1],'FontWeight','bold','FontSize',16);
set(gca,'zlim',[0,ceil(a(2))]);
colorbar('FontWeight','bold','FontSize',16);
%a=get(gca,'zlim');
%set(gca,'zlim',[0,ceil(a(2))],'clim',[0,ceil(a(2))]);

