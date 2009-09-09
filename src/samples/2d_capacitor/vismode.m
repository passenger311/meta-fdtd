function vismode(freqnumber,fignum);
% vismode(freqnumber,fignum)
% The function vismode visualises the field enhancement for different modes.
%
% Input:
% - freqnumber = integer number related to frequency in lambda2.in
% - fignum = number of figure opened, so that several modes can be visualised at the same time
%
% Files needed:
% - data.save containing the conversion factor between real world and computational units dx
% - lambda2.in containing a list of wavelengths for which the modes were calculated
% - EHT_%i.set with %i an integer number as element of the number of wavelengths
% - ../EHT_%i.set with %i same as above
%

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

filename=sprintf('EHN_%i.set',freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;
filename=sprintf('../EHN_%i.set',freqnumber);
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
tmp2 = textscan(fid, '%n %n %n %n');
tmp3 = textscan(fid1, '%n %n %n %n');
fclose(fid); fclose(fid1);

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
tmp1 = zeros(size(tmp3{1},1),size(tmp3,2));
for j = 1 : size(tmp2,2)
   tmp(:,j) = tmp2{j};
   tmp1(:,j) = tmp3{j};
end

Ez = tmp(:,1)+i*tmp(:,2);
Hz = tmp(:,3)+i*tmp(:,3);
Erefz = tmp1(:,1)+i*tmp1(:,2);
Hrefz = tmp1(:,4)+i*tmp1(:,4);

Erz = sum(Erefz)/size(Erefz,1);
Efz = reshape(Ez(:,1),floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));


x = dx*(facerange(1,1):facerange(1,3):facerange(1,2));
y = dx*(facerange(2,1):facerange(2,3):facerange(2,2));
figure(fignum);
surf(y,x,abs(Efz)/(abs(Erz)));
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

