function [E,H] = readface(freqnumber)

filename=sprintf('../dft-ref_%i.set',freqnumber);
fprintf(1,'loading file: %s\n', filename);
fid = fopen(filename, 'r');
if (fid==-1) error('file does not exist'); end;

tmp = textscan(fid, '%s',1); clear tmp;
%tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1,'CollectOutput', 1);
tmp = textscan(fid, '%n %n %n %n %n %n %n %n %n', 1);

range_tmp = zeros(1,size(tmp,2));
for j = 1 : size(tmp,2)
   range_tmp(j) = tmp{j};
end

facerange(:,:)=reshape(range_tmp,3,3)';
clear tmp, range_tmp;

tmp2 = textscan(fid, '%n %n %n %n %n %n %n %n ', 1);

tmp = zeros(1,size(tmp2,2));
for j = 1 : size(tmp2,2)
   tmp(j) = tmp2{j};
end
E = zeros(size(tmp,1),3); H=E;
E(:,1:2) = tmp(:,1:2:4)+i*tmp(:,2:2:4);
H(:,1:2) = tmp(:,5:2:8)+i*tmp(:,6:2:8);

fclose(fid);
