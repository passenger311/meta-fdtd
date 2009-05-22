function [ ] = wgsave3j( filename, e, h, range, nrange )

fid = fopen(filename,'w');
fprintf(fid,'(SET\n');

fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',nrange(1),nrange(2),1,nrange(4),nrange(5),1,nrange(7),nrange(8),1);

si = nrange(2)-nrange(1)+1;
sj = nrange(8)-nrange(7)+1;

oi = nrange(1) - range(1) ;
oj = nrange(7) - range(7) ;

for j=oj+1:oj+sj
   for i=oi+1:oi+si
      fprintf(fid,'%f %f\n',e(j,i),h(j,i));
   end
end

fprintf(fid,')SET\n');

fclose(fid);
