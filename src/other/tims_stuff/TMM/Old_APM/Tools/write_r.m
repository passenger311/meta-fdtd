function write_r(cl,ref,struc_name,prec,theta_i)

prec        =   prec-1;
fname       =   sprintf('%.0f_%02.0f_%02.0f',cl(1),cl(2),cl(3));
fname2      =   sprintf('%02.0f_%02.0f',cl(4),cl(5));

% Create Host Folder if does not Exist
FolderName = 'Outputs';
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

%Create Per-Day Folder if does not Exist
FolderName = sprintf('%s/%s',FolderName,fname);
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

% Create Structure Folder
FolderName = sprintf('%s/%s', FolderName,struc_name);
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

% Create Time Folder
FolderName = sprintf('%s/%s', FolderName,fname2);
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

FileName    =   'Inputs\input_parameters.m';

% Copy Files to Per-Day Folder if dont exist
if (exist(FileName,'file') == 0)
    copyfile(FileName, FolderName);
end

FileName    =   sprintf('Structures/%s.m',struc_name);
FileName2   =   sprintf('%s/Structure.m',FolderName);

% Copy Files to Per-Day Folder if dont exist
if (exist(FileName,'file') == 0)
    copyfile(FileName, FileName2);
end

% Save root results
FileName    =   'Reflection and transmission.txt';
FileID      =   fopen(FileName, 'w');

sym = '%';

fprintf(FileID,'Reflection and transmission Results\n\n');
a2 = sprintf('Incident angle \t%s.%ge\n\n',sym,prec);
fprintf(FileID,a2,theta_i);
fprintf(FileID,'Mode Polarisation\tReflection\tTransmission\tAbsorption\tTotal\n');
a1 = sprintf('%ss\t%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
    sym,sym,prec,sym,prec,sym,prec,sym,prec);
fprintf(FileID,a1,'TM',ref(1),ref(2),ref(3),sum(ref(1:3)));
fprintf(FileID,a1,'TE',ref(4),ref(5),ref(6),sum(ref(4:6)));
% Close Text Files
fclose(FileID);

pause(1)
% Move Files to Per-Day Folder
movefile(FileName, FolderName);

