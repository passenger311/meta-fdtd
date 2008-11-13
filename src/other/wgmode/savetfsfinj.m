function [ ] = savetfsfinj( filename, e, h, range, itrim, jtrim )

range(1) = range(1) + itrim(1);
range(2) = range(2) - itrim(2);

range(4) = range(4) + jtrim(1);
range(5) = range(5) - jtrim(2);

fid = fopen(filename,'w');
fprintf(fid,'(SET\n');

fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',range(1),range(2),1,range(4),range(5),1,range(7),range(9),1);

si = range(2)-range(1)+1;
sj = range(5)-range(4)+1;
sk = range(8)-range(7)+1;

for j=1:sj
   for i=1:si
      fprintf(fid,'%f %f\n',e(j,i),h(j,i));
   end
end

fprintf(fid,')SET\n');

fclose(fid);
