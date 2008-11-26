function saveFields (fileName,var1,var2,window_in,window_out);


fid = fopen(fileName, 'w');
fprintf(fid,'(SET\n');
fprintf(fid,'%i %i %i %i %i %i %i %i %i\n', window_out);

row_start=window_out(1)-window_in(1)+1;
row_end=window_out(2)-window_in(1)+1;
column_start=window_out(4)-window_in(4)+1;
column_end=window_out(5)-window_in(4)+1;


for column=column_start:column_end
  for row=row_start:row_end
    fprintf(fid,'%e %e\n',var1(row,column),var2(row,column));
  end
end
fprintf(fid,')SET\n');
fclose(fid);
