function sav_rt(input)
in_mat      =   input{1};
R_par_mat   =   in_mat{1};
R_per_mat   =   in_mat{2};
T_par_mat   =   in_mat{3};
T_per_mat   =   in_mat{4};
x_var       =   in_mat{5};
y_var       =   in_mat{6};
prec        =   input{2};
cl          =   input{3};
struc_name  =   input{4};
run_n       =   input{5};
run         =   input{6};


k0          =   input{7};
var1_name   =   input{8};
var2_name   =   input{9};
num_var     =   input{10};
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

FileName   =   'Inputs\input_parameters.m';

% Copy Files to Per-Day Folder if dont exist

copyfile(FileName, FolderName);


FileName    =   sprintf('Structures/%s.m',struc_name);
FileName2   =   sprintf('%s/Structure.m',FolderName);

% Copy Files to Per-Day Folder if dont exist


copyfile(FileName, FileName2);


% Create mode scan folder
FolderName = sprintf('%s/RT scan results', FolderName);
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

FileName = 'Inputs\rt_scan_in.m';

% Copy Files to mode scan Folder if dont exist

copyfile(FileName, FolderName);



% Create run folder
if run_n>1
    FolderName = sprintf('%s/Run %g', FolderName,run);
    if (exist(FolderName, 'dir') == 0)
        mkdir(FolderName);
    end
end


% Save rt results

% for loop = 1:size(R_par_mat,1)*size(R_par_mat,2)
% 
% end


FileName    =   'RT scan results.txt';
FileID      =   fopen(FileName, 'w');

sym = '%';


if num_var == 1

    a = sprintf('%s\tTM R\tTM T\tTE R\tTE T\n',var1_name);
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
    
    
    n = size(R_par_mat,2);
    
    for loop = 1:n
        par     =   x_var(loop);
        TM_R    =   real(R_par_mat(loop));
        TM_T    =   real(T_par_mat(loop));
        TE_R    =   real(R_per_mat(loop));
        TE_T    =   real(T_per_mat(loop));
        fprintf(FileID,a1,par,TM_R,TM_T,TE_R,TE_T);
        
    end

else

    a = sprintf('%s\t%s\tTM R\tTM T\tTE R\tTE T\n',var1_name,var2_name);
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);


    for loop = 1:length(x_var)
  
        for loop2 = 1:length(y_var)
            x       =   x_var(loop);
            y       =   y_var(loop2);

            TM_R    =   real(R_par_mat(loop2,loop));
            TM_T    =   real(T_par_mat(loop2,loop));
            TE_R    =   real(R_per_mat(loop2,loop));
            TE_T    =   real(T_per_mat(loop2,loop));
            fprintf(FileID,a1,x,y,TM_R,TM_T,TE_R,TE_T);
        end

    end

end

% Close Text Files
fclose(FileID);

% Copy Files to Per-Day Folder
movefile(FileName, FolderName);
