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


filename=sprintf('task_cap%i/EHN_%i.set',angnum,freqnumber);
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
tmp2 = textscan(fid, '%n %n %n %n');

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   tmp(:,j) = tmp2{j};
end

E = tmp(:,1)+i*tmp(:,2);
H = tmp(:,3)+i*tmp(:,4);

Ef = reshape(E(:,1),(floor(facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));
fclose(fid);

filename=sprintf('task_capnp%i/EHN_%i.set',angnum,freqnumber);
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
tmp2 = textscan(fid, '%n %n %n %n');

tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   tmp(:,j) = tmp2{j};
end

E = tmp(:,1)+i*tmp(:,2);
H = tmp(:,3)+i*tmp(:,4);

Ef2 = reshape(E(:,1),(floor(facerange(1,2)-facerange(1,1))/facerange(1,3)+1),floor((facerange(2,2)-facerange(2,1))/facerange(2,3)+1));
fclose(fid);

x = 10*(facerange(1,1):facerange(1,3):facerange(1,2));
y = 10*(facerange(2,1):facerange(2,3):facerange(2,2));
maxim = 1.05*max(max([abs(Ef),abs(Ef2)]));

figure(fignum);
hold on;
surf(y,x,abs(Ef));
line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('y in nm','FontWeight','bold','FontSize',30);
ylabel('x in nm','FontWeight','bold','FontSize',30);
title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',30);
shading interp;
view(90,-90);
%axis equal;
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',30);
axis([-Inf Inf -Inf Inf 0 maxim 0 maxim]);
colorbar('FontWeight','bold','FontSize',30);
hold off;
filename = sprintf('wonp_%i_%i.png',angnum,lambda(freqnumber));
print(gcf, '-dpng', '-r300', filename);

figure(fignum+1);
hold on;
surf(y,x,abs(Ef2));
line([-190,-190,-210,-210,-190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
line([190,190,210,210,190],[-1500,1500,1500,-1500,-1500],'Color',[1 1 1],'Linewidth',1);
xlabel('y in nm','FontWeight','bold','FontSize',30);
ylabel('x in nm','FontWeight','bold','FontSize',30);
title(sprintf('Wavelength %i nm',lambda(freqnumber)),'FontWeight','bold','FontSize',30);
shading interp;
view(90,-90);
foo = get(gca,'dataaspectratio');
set(gca,'dataaspectratio',[foo(1) foo(1) foo(3)],'FontWeight','bold','FontSize',30);
axis([-Inf Inf -Inf Inf 0 maxim 0 maxim]);
colorbar('FontWeight','bold','FontSize',30);
%hold off;
filename = sprintf('wnp_%i_%i.png',angnum,lambda(freqnumber));
print(gcf, '-dpng', '-r300', filename);

