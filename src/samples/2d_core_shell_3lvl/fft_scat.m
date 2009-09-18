s = {'x','y','z'};

fid = fopen('data.save','r');
tline = fgetl(fid); 
dx = sscanf(tline,'%e');
tline = fgetl(fid);
boxsize = sscanf(tline,'%e %e %e');
fclose(fid);
filename=sprintf('s_scat.pspec');
fid3 = fopen(filename,'w');
if (fid3==-1) continue; end;
fprintf(1,'file: %s\n', filename);
for i = 1:2
  filename=sprintf('fft+%sscat.pspec',s{i});
  fid = fopen(filename, 'r');
  filename=sprintf('fft-%sscat.pspec',s{i});
  fid2 = fopen(filename,'r');
  if (fid2==-1) continue; end;
  for j=1:7
      tline = fgetl(fid);
      tline = fgetl(fid2);
      if i==1
        fprintf(fid3,tline);
        fprintf(fid3,'\n');
      end
  end
  tmp = zeros(3,3); wg=zeros(3,3);np=zeros(3,3);
  k=1;
  while feof(fid) == 0
      tline = fgetl(fid);
      wg(i,:) = sscanf(tline, '%e %e %e');
      tline = fgetl(fid2);
      np(i,:) = sscanf(tline, '%e %e %e');
      tmp(i,1) = 0;
      lambda(i,k) = np(i,1);
      for j = 2:1:3; tmp(i,j) = wg(i,j)-np(i,j); end;
      stot(i,:,k) = tmp(i,:)*(boxsize(mod(i,3)+1))*(boxsize(mod(i+1,3)+1));
      k=k+1;
  end;
  k_max=k;
  fclose(fid);
  fclose(fid2);
end
stotal = sum(sum(stot,1),2);
for k=1:k_max-1
  fprintf(fid3,'%e %e\n', dx/lambda(2,k),dx*stotal(1,1,k));
end
fclose(fid3);
filename=sprintf('../fft_ref-zabs.pspec');
fid = fopen(filename, 'r');
for j=1:7
  tline = fgetl(fid);
end
sz=zeros(1,3);
k=1;
while feof(fid) == 0
  tline = fgetl(fid);
  sz(1,:) = sscanf(tline, '%e %e %e');
  s_inj(1,1,k) = sz(1,2);
  k=k+1;
end
fclose(fid);
np_abs = stotal ./ s_inj;
filename=sprintf('s_relscat.pspec');
fid3 = fopen(filename,'w');
for k=1:k_max-1
  fprintf(fid3,'%e %e\n', dx/lambda(2,k),dx*abs(np_abs(k)));
end
fclose(fid3);
exit
