function [ ] = wgsave3k( filename, e1, e2, e3, range, nrange )

fid = fopen(filename,'w');
fprintf(fid,'(SET\n');

fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',nrange(1),nrange(2),1,nrange(4),nrange(5),1,nrange(7),nrange(8),1);

si = nrange(2)-nrange(1)+1;
sj = nrange(5)-nrange(4)+1;

oi = nrange(1) - range(1) ;
oj = nrange(4) - range(4) ;

for j=oj+1:oj+sj
   for i=oi+1:oi+si
      fprintf(fid,'%f %f %f\n',e1(j,i),e2(j,i),e3(j,i));
   end
end

fprintf(fid,')SET\n');

fclose(fid);
