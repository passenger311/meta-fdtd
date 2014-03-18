% function roughness2
close all
clear all
clc

% Get geo file names

list    =   dir('geo_*.in');
n       =   0;
num     =   length(list);
nm      =   cell(2,num);
fnm     =   cell(1,num);
number  =   fnm;

for n   =   1:num
    
    a       =   list(n).name;
    st      =   find(a == '_');
    fi      =   find(a == '.');
    a       =   a(st(1)+1:fi-1);
    nm{1,n}     =   a;
    
end


% Compare geo_filling names with geo names and assign the same 
% random seed to the same interfaces.

a   =   0;
for n   =   1:num
    
    ind     =   find(nm{1,n} == 'f', 1);
    
    if isempty(ind) == 0
        a       =   a+1;
        ind2    =   find(nm{1,n} == '_');
        fnm     =   nm{1,n}(ind2(1)+1:ind2(2)-1);
        nm1     =   nm(1,:);
        ind3    =   strfind(nm1,fnm);
        ind4    =   find(not(cellfun('isempty',ind3)));
        nm{2,ind4(1)}   =   a;
        nm{2,ind4(2)}   =   a;
    end
    
end

ind     =   find(cellfun('isempty',nm));

for n   =   1:length(ind)
    nm{ind(n)}      =   ind(n);
end


% Read interface file

filename    =   'Interfaces.txt';
fileID      =   fopen(filename,'r');
loop        =   0;

while 1
    
    tline   =   fgetl(fileID);
    
    if ischar(tline) == 0
        break
    end
    loop    =   loop+1;
    inf     =   textscan(tline,'%s');
    name    =   inf{1}{1};
    layer.(name)    =   loop;
    number{loop}    =   name;
    pos.(name)  =   str2double(inf{1}{2});
    rpos.(name) =   str2double(inf{1}{3});
    eps.(name)  =   str2double(inf{1}{4});
    filling.(name)  =   str2double(inf{1}{5});
    
end



for n = 1:num
    
    filename    =   sprintf('geo_%s.in',nm{1,n});
    fileID      =   fopen(filename,'r');
    
    for loop = 1:4
        tline   =   fgetl(fileID);
    end
    
    fclose(fileID);
    
    ln      =   nm{1,n};
    ind     =   find(ln == '_');
    
    if isempty(ind) == 0
        ln  =   ln(ind(1)+1:ind(2)-1);
    end
 
    sz      =   textscan(tline,'%f');
    xle     =   (sz{1}(2)-sz{1}(1))+1;  % Length of interface
    yle     =   (sz{1}(5)-sz{1}(4))+1;  % Height of interface
    yst     =   sz{1}(4);
    
    lind    =   layer.(ln);
    low     =   number{lind+1};
    
    if isempty(ind) == 0
        eps1    =   filling.(ln);
        eps2    =   filling.(low);
    else
        eps1    =   eps.(ln);
        eps2    =   eps.(low);  
    end
    
    % Create interface
    
    intfc   =   ones(1,(xle*2)-1)*pos.(ln);
    x       =   1:0.5:xle;
    
    out_mat =   cell(1,yle*2);
    
    for loop    =   1:yle*2
        
        int     =   intfc-(0.5*(loop-1)+yst);
        si      =   sign(int)+1;
        diff    =   si(2:end)-si(1:end-1);
        ind     =   find(abs(diff)==2);
        x_ex    =   zeros(1,length(ind));
        
        for var = 1:length(ind)
            y1  =   int(ind(var));
            y2  =   int(ind(var)+1);
            x1  =   x(ind(var));
            x2  =   x(ind(var)+1);
            m   =   (y2-y1)/(x2-x1);
            c   =   y2-(m*x2);
            z_p =   -c/m;
            x_ex(var) = z_p;
        end
        
        x_mat_zero  =   [x x_ex];
        int_zero    =   [int zeros(1,var)];
        
        int2    =   int-0.5;
        si      =   sign(int2)+1;
        diff    =   si(2:end)-si(1:end-1);
        ind     =   find(abs(diff)==2);
        x_ex    =   zeros(1,length(ind));
        
        for var = 1:length(ind)
            y1  =   int2(ind(var));
            y2  =   int2(ind(var)+1);
            x1  =   x(ind(var));
            x2  =   x(ind(var)+1);
            m   =   (y2-y1)/(x2-x1);
            c   =   y2-(m*x2);
            z_p =   -c/m;
            x_ex(var) = z_p;
        end
        
        x_mat_top   =   [x_mat_zero x_ex];
        int_top     =   [int_zero ones(1,var)*0.5];
        
        out     =   [x_mat_top;int_top];
        out     =   out';
        s_out   =   sortrows(out,1);
        
        ind1    =   find(s_out(:,2)>=0.5);
        ind2    =   find(s_out(:,2)<=0);
        s_out(ind1,2)   =   0.5;
        s_out(ind2,2)   =   0;
        
        out_mat{loop}   =   s_out;
       
    end
    
    
    
end

figure(1)
pcolor(int)
shading flat
