function out = read_custom

FileName    =   'custom.txt';

FileID      =   fopen(FileName, 'r');
loop = 0;

    while 1
        loop = loop+1;
        tline = fgetl(FileID);
        if ~ischar(tline),   break, end
        
        if isempty(tline) == 0
            custom_mat(loop) = str2num(tline)*1e9;
        end
    end

    out = custom_mat;
    
   