function out = sav_mode(out_mat,prec,cl,struc_name,run_n,run,num,mode,k0,var1_name,var2_name)

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
FolderName = sprintf('%s/Mode scan results', FolderName);
if (exist(FolderName, 'dir') == 0)
    mkdir(FolderName);
end

FileName = 'Inputs\m_scan_in.m';

% Copy Files to mode scan Folder if dont exist

    copyfile(FileName, FolderName);



% Create run folder
if run_n>1
    FolderName = sprintf('%s/Run %g', FolderName,run);
    if (exist(FolderName, 'dir') == 0)
        mkdir(FolderName);
    end
end

% Create mode folder
if num>1
    FolderName = sprintf('%s/Mode %g', FolderName,mode);
    if (exist(FolderName, 'dir') == 0)
        mkdir(FolderName);
    end
end

% Save root results
FileName    =   'Mode scan results.txt';
FileID      =   fopen(FileName, 'w');

sym = '%';
 fprintf(FileID,'Modal Results\n\n');
  a = sprintf('k0 \t%s.%ge\n\n',sym,prec);
 fprintf(FileID,a,k0);
if size(out_mat,2) == 6
    a = sprintf('%s\tRe(z)\tIm(z)\tRe(R)\tIm(R)\tRe(S)\tIm(S)\tk0\tVg\n',var1_name);
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
    n = size(out_mat,1);
    for loop = 1:n
        par     =   out_mat(loop,1);
        z_re    =   real(out_mat(loop,2));
        z_im    =   imag(out_mat(loop,2));
        R_re    =   real(out_mat(loop,3));
        R_im    =   imag(out_mat(loop,3));        
        S_re    =   real(out_mat(loop,4));
        S_im    =   imag(out_mat(loop,4));
        k0      =   out_mat(loop,5);
        vg      =   out_mat(loop,6);
        fprintf(FileID,a1,par,z_re,z_im,R_re,R_im,S_re,S_im,k0,vg);
        
    end
    
else
    a = sprintf('%s\t%s\tRe(z)\tIm(z)\tRe(R)\tIm(R)\tRe(S)\tIm(S)\tk0\tVg\n',var1_name,var2_name);
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
    n = size(out_mat,1);
    for loop = 1:n
        par1    =   out_mat(loop,1);
        par2    =   out_mat(loop,2);
        z_re    =   real(out_mat(loop,3));
        z_im    =   imag(out_mat(loop,3));
        R_re    =   real(out_mat(loop,4));
        R_im    =   imag(out_mat(loop,4));        
        S_re    =   real(out_mat(loop,5));
        S_im    =   imag(out_mat(loop,5));
        k0      =   out_mat(loop,6);
        vg      =   out_mat(loop,7);
        fprintf(FileID,a1,par1,par2,z_re,z_im,R_re,R_im,S_re,S_im,k0,vg);
        
    end
    
end

% Close Text Files
fclose('all');

% Copy Files to Per-Day Folder

movefile(FileName, FolderName,'f');

% Save root results
FileName    =   'Mode scan outputs.txt';
FileID      =   fopen(FileName, 'w');
c   =   299792458;
lambdap =   (2*pi*c);
sym = '%';
fprintf(FileID,'Modal Results\n\n');
a = sprintf('k0 \t%s.%ge\n\n',sym,prec);
%fprintf(FileID,a,k0);
if size(out_mat,2) == 6
    a = sprintf('%s\tRe(Neff)\tIm(Neff)\tAbsorp\tRe(k)\tIm(k)\tBound\tVg\n',var1_name);
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
    n = size(out_mat,1);
    
    for loop = 1:n
        par     =   out_mat(loop,1);%.*1.62./11.22e6./2.9989e8;%/lambdap;
        z       =   out_mat(loop,2);
        R       =   out_mat(loop,3);
        S       =   out_mat(loop,4);
        k0      =   out_mat(loop,5);
        vg      =   (out_mat(loop,6));%-0.01);
        w       =   exp(2*z)+((R^2)*exp(-2*z))+S;
        n       =   sqrt(w);
        n       =   n*sign(imag(n));
        if loop == 1
            oldn = n;
        end
        if sign(real(n)/real(oldn)) == -1
            d = real(n)-real(oldn);
            d = abs(round(d/real(n)));
            if d == 2
                n = -n;
            end
        end
        k       =   n*k0;%./11.22e6.*1.62;
        absorp  =   (0.02*k0*imag(n));
        gammac  =   k0.*(exp(z)+R.*exp(-z));
        gammas  =   k0.*(exp(z)-R.*exp(-z));
        if real(gammac) > 0 && real(gammas) > 0
            b = 1;
        else
            b = 0;
        end
%         prop    =   1/(2*k0*abs(imag(n)));
%         sdc     =   1/(2*abs(imag(gammac)));
%         sds     =   1/(2*abs(imag(gammas)));
        fprintf(FileID,a1,par,real(n),imag(n),absorp,real(k),imag(k),b,vg);
        %fprintf(FileID,a1,par,real(n),imag(n),absorp,prop,b,sdc,sds,vg);
        k_mat(loop) = k;
        par_mat(loop) = par;
        oldn = n;
        
    end
    
else
    a = sprintf('%s\t%s\tRe(z)\tIm(z)\tRe(R)\tIm(R)\tRe(S)\tIm(S)\tVg\n',var1_name,var2_name);
    %fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec,sym,prec);
    n = size(out_mat,1);
    for loop = 1:n
        par1    =   out_mat(loop,1);
        par2    =   out_mat(loop,2);
        z       =   out_mat(loop,3);
        R       =   out_mat(loop,4);
        S       =   out_mat(loop,5);
        k0      =   out_mat(loop,6);
        vg      =   out_mat(loop,7);
        w       =   exp(2*z)+((R^2)*exp(-2*z))+S;
        n       =   sqrt(w);
        n       =   n*sign(imag(n));
        absorp  =   0.02*k0*imag(n);
        gammac  =   k0.*(exp(z)+R.*exp(-z));
        gammas  =   k0.*(exp(z)-R.*exp(-z));
        if real(gammac) > 0 && real(gammas) > 0
            b = 1;
        else
            b = 0;
        end
        prop    =   1/(2*k0*abs(imag(n)));
        sdc     =   1/(2*abs(imag(gammac)));
        sds     =   1/(2*abs(imag(gammas)));
        
        fprintf(FileID,a1,par1,par2,real(n),imag(n),absorp,prop,b,sdc,sds,vg);

        
    end
    
end

% Close Text Files
fclose(FileID);

% Copy Files to Per-Day Folder
movefile(FileName, FolderName);
out = [k_mat; par_mat];
