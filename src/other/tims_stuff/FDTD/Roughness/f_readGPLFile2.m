function output = f_readGPLFile2( fileName )
    
    fileID = fopen(fileName,'r');
    %fprintf(1,'\nReading from data file...');
    for i = 1:13
        tline{i} = fgetl(fileID);
    end
    txt=sprintf('%s\t',tline{:});
    params = str2num(tline{11}(4:end));
    rows = length((params(1):params(3):params(2)));
    cols = 1;
    skip = length(txt);
    fclose(fileID);
    fileID = fopen(fileName,'r');
    freadData = fread(fileID,'char=>char');
    freadData = sscanf(freadData(skip:end),'%g');
    l = length(freadData);
    rows = l/cols;
    freadData = reshape(freadData,cols,rows)';
    fclose(fileID);
    params(2) = params(3)*rows-1;


    energy = freadData;
   % fprintf(1,'\nComplete\n');
    output = {params,energy};
    
end
