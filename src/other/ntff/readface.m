function [facerange,E,H] = readface(freqnumber,face)

%open file containing E and H fields on face
a = {'+x','-x','+y','-y','+z','-z'};
filename=sprintf('dft%s_%i.set',a{face},freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;

%read range of box
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);
range_tmp = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2)
   range_tmp(j) = tmp{j};
end
facerange(:,:)=reshape(range_tmp,3,3)'; clear tmp;


%read complex E and H fields
tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n');
tmp = zeros(size(tmp2{1},1),size(tmp2,2));
for j = 1 : size(tmp2,2)
   tmp(:,j) = tmp2{j};
end
E = zeros(size(tmp,1),3); H=E;
if ((face ==1) | (face == 2)) 
   E(:,2:3) = tmp(:,1:2:4)+i*tmp(:,2:2:4);
   H(:,2:3) = tmp(:,5:2:8)+i*tmp(:,6:2:8);
elseif ((face ==3) | (face == 4))
   E(:,1) = tmp(:,1)+i*tmp(:,2);
   E(:,3) = tmp(:,3)+i*tmp(:,4);
   H(:,1) = tmp(:,5)+i*tmp(:,6);
   H(:,3) = tmp(:,7)+i*tmp(:,8);
elseif ((face ==5) | (face == 6))
   E(:,1:2) = tmp(:,1:2:4)+i*tmp(:,2:2:4);
   H(:,1:2) = tmp(:,5:2:8)+i*tmp(:,6:2:8);
end

fclose(fid);
