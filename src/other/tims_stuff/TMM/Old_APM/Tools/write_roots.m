function write_roots(cl,roots,struc_name,prec,R,S,k0)

prec        =   prec-1;
fname       =   sprintf('%.0f_%02.0f_%02.0f',cl(1),cl(2),cl(3));
fname2      =   sprintf('%02.0f_%02.0f',cl(4),cl(5));
num         =   size(roots,1);
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

% Copy Files to Per-Day Folder
copyfile(FileName, FolderName);

FileName    =   sprintf('Structures/%s.m',struc_name);
FileName2   =   sprintf('%s/Structure.m',FolderName);

% Copy Files to Per-Day Folder
copyfile(FileName, FileName2);


% Save root results
FileName    =   'Modes.txt';
FileID      =   fopen(FileName, 'w');

n = size(roots,1);
sym = '%';

fprintf(FileID,'Modal Results\n\n');
a2 = sprintf('k0 \t%s.%ge\nR \t%s.%ge\nS \t%s.%ge\nNumber of roots \t%s.%gg\n\n',sym,prec,sym,prec,sym,prec,sym,prec);
fprintf(FileID,a2,k0,R,S,num);
fprintf(FileID,'Root\tRe(z)\tIm(z)\tRe(Neff)\tIm(Neff)\tAbsorp(1/cm)\tBound\n');
a1 = sprintf('%s.%gg\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%gg\n',...
    sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);

for loop = 1:n
    fprintf(FileID,a1,loop,real(roots(loop,1)),...
        imag(roots(loop,1)),roots(loop,2),roots(loop,3),roots(loop,4),roots(loop,5));
end

% Close Text Files
fclose(FileID);

% Copy Files to Per-Day Folder
copyfile(FileName, FolderName);

% Move to tmp Folder
FolderName = 'tmp';
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end
movefile(FileName,FolderName)
FileName    =   'Mode scan results.txt';
FileID      =   fopen(FileName, 'w');

sym = '%';
% fprintf(FileID,'Modal Results\n\n');
% a = sprintf('k0 \t%s.%ge\n\n',sym,prec);
% fprintf(FileID,a,k0);

    a = sprintf('k0\tRe(z)\tIm(z)\tRe(R)\tIm(R)\tRe(S)\tIm(S)\tk0\tVg\n');
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
   
    for loop = 1:n
        
        z_re    =   real(roots(loop,1));
        z_im    =   imag(roots(loop,1));
        R_re    =   real(R);
        R_im    =   imag(R);        
        S_re    =   real(S);
        S_im    =   imag(S);
        
        vg      =   1;
        fprintf(FileID,a1,k0,z_re,z_im,R_re,R_im,S_re,S_im,k0,vg);
        
    end



