function vismodecomp(angnum,freqnumber,fignum);

fid1 = fopen('lambda2.in', 'r');
if (fid1==-1) error('file lambda.in does not exits'); end;
k=1;
while feof(fid1) ==0
   tline = fgetl(fid1);
   lambda(k)=sscanf(tline,'%e');
   k=k+1;
end
fclose(fid1);


filename=sprintf('task_cap%i/F_%i.set',angnum,freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;

%tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1,'CollectOutput', 1);
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2)
   range(j) = tmp{j};
end
facerange(:,:)=reshape(range,3,3)';
clear tmp;
%tmp = textscan(fid,'%n %n %n %n %n %n %n %n %n %n %n %n', 'CollectOutput', 1);
%tmp = tmp{1};
tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n %n ');

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
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

Ef = reshape(Etmp,floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));
Hf = reshape(Htmp,floor((facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));

x = 2*(facerange(1,1):facerange(1,3):facerange(1,2));
y = 2*(facerange(2,1):facerange(2,3):facerange(2,2));

fclose(fid);

filename=sprintf('F_%i.set',freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;

%tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1,'CollectOutput', 1);
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2)
   range(j) = tmp{j};
end
facerange(:,:)=reshape(range,3,3)';
clear tmp;
%tmp = textscan(fid,'%n %n %n %n %n %n %n %n %n %n %n %n', 'CollectOutput', 1);
%tmp = tmp{1};
tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n %n %n %n %n ');

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
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

Eref = sum(sqrt(E(:,1).*conj(E(:,1))+E(:,2).*conj(E(:,2))+E(:,3).*conj(E(:,3))))/size(E(:,1),1);
Href = sum(sqrt(H(:,1).*conj(H(:,1))+H(:,2).*conj(H(:,2))+H(:,3).*conj(H(:,3))))/size(H(:,1),1);

fclose(fid);

Ef = Ef./Eref;
Hf = Hf./Href;

maxim = 1.05*max(max(Ef));
maximH = 1.05*max(max(Hf));

figure(fignum);
clf reset
hold on;
surf(y,x,Ef);
%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('y in nm','FontWeight','bold','FontSize',20);
ylabel('x in nm','FontWeight','bold','FontSize',20);
title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',24);
shading interp;
view(90,-90);
%axis equal;
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',20);
axis([-Inf Inf -Inf Inf 0 maxim 0 maxim]);
colorbar('South','FontWeight','bold','FontSize',20);
hold off;
filename = sprintf('E_%ideg_%inm.png',5*(angnum-1),lambda(freqnumber));
print(gcf, '-dpng', '-r200', filename);

figure(fignum+1);
clf reset
hold on;
surf(y,x,Hf);
%line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
%line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('y in nm','FontWeight','bold','FontSize',20);
ylabel('x in nm','FontWeight','bold','FontSize',20);
title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',24);
shading interp;
view(90,-90);
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',20);
axis([-Inf Inf -Inf Inf 0 maximH 0 maximH]);
colorbar('South','FontWeight','bold','FontSize',20);
hold off;
filename = sprintf('H_%ideg_%inm.png',5*(angnum-1),lambda(freqnumber));
print(gcf, '-dpng', '-r200', filename);

