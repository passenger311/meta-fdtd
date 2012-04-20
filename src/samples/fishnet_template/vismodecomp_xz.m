function vismodecomp_xz(file,freqnumber,fignum);

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



filename=sprintf('F_%s_%i.set',file,freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;

%tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1,'CollectOutput', 1);
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2);
   range(j) = tmp{j};
end
facerange(:,:)=reshape(range,3,3)';
clear tmp;
%tmp = textscan(fid,'%n %n %n %n %n %n %n %n %n %n %n %n', 'CollectOutput', 1);
%tmp = tmp{1};
tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n %n ');

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

Etmp = sqrt(E(:,1).*conj(E(:,1))+E(:,2).*conj(E(:,2))+E(:,3).*conj(E(:,3)));
Htmp = sqrt(H(:,1).*conj(H(:,1))+H(:,2).*conj(H(:,2))+H(:,3).*conj(H(:,3)));

Ef = reshape(E,floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(3,2)-facerange(3,1))/facerange(3,3)+1),3);
%Ef = reshape(Etmp,floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));
%Hf = reshape(Htmp,floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));

Ef = Ef./max(max(max(abs(Ef))));

%size(Ef)
%size(Hf)
angle(max(max(Ef(:,:,1))));
[val1,pos1]=max(abs(Ef(:,:,1)));
[val2,pos2]=max(val1);
angle(Ef(pos1(pos2),pos2,1));
%Efx = real(Ef.*exp(-i*angle(max(max(Ef(:,:,1))))));%2.7349));%1.532));%-2.0718));%.*exp(-i*angle(max(max(Ef(:,:,1))))));
Efx = real(Ef.*exp(-i*angle(max(max(Ef(:,:,1))))));
Efz = real(Ef.*exp(-i*angle(max(max(Ef(:,:,3))))));
Ef = sqrt(Ef.*conj(Ef));

%Ef = imag(Ef);
x = (facerange(1,1):facerange(1,3):facerange(1,2)).*dx;
y = (facerange(2,1):facerange(2,3):facerange(2,2)).*dx;
z = (facerange(3,1):facerange(3,3):facerange(3,2)).*dx;

fclose(fid);

maxim = ceil(max(max(max(Ef))));
%maxim = 1.1*maxim;
%maxim=40
%maximH = 1.05*max(max(Hf));

Ef=[Ef; Ef(1,:,:)];
Efx=[Efx; Efx(1,:,:)];
Efz=[Efz; Efz(1,:,:)];
a = size(Ef,1);
b = size(Ef,2);
x=[x, -x(1)];

figure(fignum);
clf reset
hold on;
E_plot = sqrt(Ef(1:a,1:b,1).*Ef(1:a,1:b,1)+Ef(1:a,1:b,2).*Ef(1:a,1:b,2)+Ef(1:a,1:b,3).*Ef(1:a,1:b,3));
surf(z(1:b),x(1:a),E_plot);%sqrt(Ef(2:a,2:b,1).*conj(Ef(2:a,2:b,1))));
%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('z in nm','FontWeight','bold','FontSize',16);
ylabel('x in nm','FontWeight','bold','FontSize',16);
title(sprintf('Wavelength %4.1f nm',lambda(freqnumber)),'FontWeight','bold','FontSize',16);
shading interp;
view(90,-90);
%axis equal;
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',16);
axis([-Inf Inf -Inf Inf -1 maxim+1 0 maxim]);
%caxis([0 maxim])
colorbar('EastOutside','FontWeight','bold','FontSize',16);
hold off;
filename = sprintf('F_%s_%4.1fnm.png',file,lambda(freqnumber));
print(gcf, '-dpng', '-r100', filename);
%filename = sprintf('F_%s_%4.1fnm.eps',file,lambda(freqnumber));
%print(gcf, '-depsc', '-r300', filename);
close(gcf);

maxim = max(max(max(abs(Efx(:,:,1)))));
%maxim = 1.1*maxim;
figure(fignum+1);
clf reset
hold on;
E_plot = Efx(1:a,1:b,1);
surf(z(1:b),x(1:a),E_plot);%sqrt(Ef(2:a,2:b,1).*conj(Ef(2:a,2:b,1))));
%%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('z in nm','FontWeight','bold','FontSize',16);
ylabel('x in nm','FontWeight','bold','FontSize',16);
title(sprintf('Wavelength %4.1f nm',lambda(freqnumber)),'FontWeight','bold','FontSize',16);
shading interp;
view(90,-90);
%%axis equal;
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',16);
axis([-Inf Inf -Inf Inf -maxim-1 maxim+1 -maxim maxim]);
colorbar('EastOutside','FontWeight','bold','FontSize',16);
hold off;
filename = sprintf('Fx_%s_%4.1fnm.png',file,lambda(freqnumber));
print(gcf, '-dpng', '-r100', filename);
close(gcf)


%maxim = max(max(max(abs(Efx(:,:,2)))));
%maxim = 1.1*maxim;
%figure(fignum+2);
%clf reset
%hold on;
%E_plot = Efx(2:a,2:b,2);
%surf(y(2:b),x(2:a),E_plot);%sqrt(Ef(2:a,2:b,1).*conj(Ef(2:a,2:b,1))));
%%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%xlabel('y in nm','FontWeight','bold','FontSize',16);
%ylabel('x in nm','FontWeight','bold','FontSize',16);
%title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',16);
%shading interp;
%view(90,-90);
%%axis equal;
%foo = get(gca,'dataaspectratio');
%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',16);
%axis([-Inf Inf -Inf Inf -maxim maxim -maxim maxim]);
%colorbar('EastOutside','FontWeight','bold','FontSize',16);
%hold off;

maxim = max(max(max(abs(Efz(:,:,3)))));
%maxim = 1.1*maxim;
figure(fignum+2);
clf reset
hold on;
E_plot = Efz(1:a,1:b,3);
surf(z(1:b),x(1:a),E_plot);%sqrt(Ef(2:a,2:b,1).*conj(Ef(2:a,2:b,1))));
%%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('z in nm','FontWeight','bold','FontSize',16);
ylabel('x in nm','FontWeight','bold','FontSize',16);
title(sprintf('Wavelength %4.1f nm',lambda(freqnumber)),'FontWeight','bold','FontSize',16);
shading interp;
view(90,-90);
%%axis equal;
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',16);
axis([-Inf Inf -Inf Inf -maxim-1 maxim+1 -maxim maxim]);
colorbar('EastOutside','FontWeight','bold','FontSize',16);
hold off;
filename = sprintf('Fz_%s_%4.1fnm.png',file,lambda(freqnumber));
print(gcf, '-dpng', '-r100', filename);
close(gcf)

%%filename = sprintf('E_%s_task%i_%inm.eps',file,angnum,lambda(freqnumber));
%%print(gcf, '-depsc', '-r300', filename);

%figure(fignum+1);
%clf reset
%hold on;
%surf(y,x,Hf);
%%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%xlabel('y in nm','FontWeight','bold','FontSize',16);
%ylabel('x in nm','FontWeight','bold','FontSize',16);
%title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',16);
%shading interp;
%view(90,-90);
%foo = get(gca,'dataaspectratio');
%set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',16);
%axis([-Inf Inf -Inf Inf 0 maximH 0 maximH]);
%colorbar('East','FontWeight','bold','FontSize',16);
%hold off;
%filename = sprintf('H_%ideg_%inm.png',5*(angnum-1),lambda(freqnumber));
%print(gcf, '-dpng', '-r200', filename);




filename=sprintf('Fx_%s_%i.vtk',file,freqnumber);
fprintf(1,'writing file: %s\n', filename);
fid = fopen(filename, 'w');

fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'META: 0.1.100323, 3D\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET RECTILINEAR_GRID\n');
fprintf(fid,'DIMENSIONS %i %i %i\n',size(x,2),1,size(z,2));
fprintf(fid,'X_COORDINATES %i float\n',size(x,2));
fprintf(fid,'  %f',x);
fprintf(fid,'\nY_COORDINATES %i float\n',1);
fprintf(fid,'  %f',y);
fprintf(fid,'\nZ_COORDINATES %i float\n',size(z,2));
fprintf(fid,'  %f',z);
fprintf(fid,'\nPOINT_DATA %i\n', size(x,2)*size(z,2));
fprintf(fid,'VECTORS vectors float\n');
for j = 1:b;
  for k = 1:a;
    fprintf(fid,'%e %e %e\n',Efx(k,j,1),Efx(k,j,2),Efx(k,j,3));
  end
end
fclose(fid);

filename=sprintf('Fz_%s_%i.vtk',file,freqnumber);
fprintf(1,'writing file: %s\n', filename);
fid = fopen(filename, 'w');

fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'META: 0.1.100323, 3D\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET RECTILINEAR_GRID\n');
fprintf(fid,'DIMENSIONS %i %i %i\n',size(x,2),1,size(z,2));
fprintf(fid,'X_COORDINATES %i float\n',size(x,2));
fprintf(fid,'  %f',x);
fprintf(fid,'\nY_COORDINATES %i float\n',size(y,2));
fprintf(fid,'  %f',y);
fprintf(fid,'\nZ_COORDINATES %i float\n',size(z,2));
fprintf(fid,'  %f',z);
fprintf(fid,'\nPOINT_DATA %i\n', size(x,2)*size(z,2));
fprintf(fid,'VECTORS vectors float\n');
for j = 1:b
  for k = 1:a
    fprintf(fid,'%e %e %e\n',Efz(k,j,1),Efz(k,j,2),Efz(k,j,3));
  end
end
fclose(fid);
