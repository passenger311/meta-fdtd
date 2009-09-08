fid = fopen('data.save','r');
tline = fgetl(fid);
dx = sscanf(tline, '%e');
fclose(fid);

filename=sprintf('../fft_ref-zabs.pspec');
fid = fopen(filename, 'r');
filename=sprintf('fft-zscat.pspec');
fid2 = fopen(filename,'r');
filename=sprintf('fft+zabs.pspec');
fid3 = fopen(filename,'r');
filename=sprintf('fft_R.pspec');
fid4 = fopen(filename,'w');
filename=sprintf('fft_T.pspec');
fid5 = fopen(filename,'w');
filename=sprintf('fft_A.pspec');
fid6 = fopen(filename,'w');
fprintf(1,'file: %s\n', filename);
for j=1:7
   tline = fgetl(fid);
   tline = fgetl(fid2);
   tline = fgetl(fid3);      
   fprintf(fid4,tline);
   fprintf(fid4,'\n');
   fprintf(fid5,tline);
   fprintf(fid5,'\n');
   fprintf(fid6,tline);
   fprintf(fid6,'\n');
end
tmp = zeros(1,3); tmp2 = zeros(1,3); yscat = zeros(1,3); yabs = zeros(1,3);
while feof(fid) == 0
   tline = fgetl(fid);
   ref = sscanf(tline, '%e %e %e');
   tline = fgetl(fid2);
   yscat = sscanf(tline, '%e %e %e');
   tline = fgetl(fid3);
   yabs = sscanf(tline, '%e %e %e');
   tmp(1) = dx/yscat(1); tmp2(1) = tmp(1);
   for j = 2:1:3; tmp(j) = -yscat(j)/ref(j); end;
   fprintf(fid4,'%e %e %e\n', tmp);
   for j = 2:1:3; tmp2(j) = yabs(j)/ref(j); end;
   fprintf(fid5,'%e %e %e\n', tmp2);
   tmp2(1)=0;tmp(1)=1-tmp(1);
   fprintf(fid6,'%e %e %e\n', 1-(tmp2+tmp));
end;
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
exit
