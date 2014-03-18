close all
clear all
clc

%cd 28132.localhost/Data/
layername   =   'layer1';
filename    =   sprintf('geo_%s.in',layername);
fileID      =   fopen(filename,'r');


loop    =   0;

while 1
    
    loop = loop+1;
    
    if loop == 5
        sz(1,:)    =   str2num(tline);
    elseif loop > 4
        fields(loop-5,:) = str2num(tline);
    end
    
    tline = fgetl(fileID);
    
    if isequal(tline,')SET'), break,end
    
    
end
fclose(fileID);
field_nam = {'Ex','Ey','Ez','Hx','Hy','Hz'};



xle = (sz(2)-sz(1))+1;
yle = (sz(5)-sz(4))+1;
y   =   1:yle;
x   =   1:xle;

% for loop = 1:6
% f = fields(:,loop);
% fr = reshape(f,xle,yle);
% fi.(field_nam{loop}) = fr;
% a = field_nam{loop}
% fiel.(a) = fr;
% fr = fr';
% fr = flipud(fr);
% fr(:,1:3)
% %figure(loop)
% %pcolor(x,y,fr')
% %shading flat
% %title(field_nam{loop})
% end
% fr1     =   fr';
% fr_r = reshape(fr,xle*yle,1);
% 
% isequal(fr_r,f)

filename    =   sprintf('geo_%s_new.in',layername);
fileID      =   fopen(filename,'r');


loop    =   0;

while 1
    
    loop = loop+1;
    
    if loop == 5
        sz(1,:)    =   str2num(tline);
    elseif loop > 4
        fields2(loop-5,:) = str2num(tline);
    end
    
    tline = fgetl(fileID);
    
    if isequal(tline,')SET'), break,end
    
    
end
fclose(fileID);
field_nam = {'Ex','Ey','Ez','Hx','Hy','Hz'};



xle = (sz(2)-sz(1))+1;
yle = (sz(5)-sz(4))+1;
y   =   1:yle;
x   =   1:xle;

for loop = 1:6
f = fields2(:,loop);
fr = reshape(f,xle,yle);
fi.(field_nam{loop}) = fr;
a = field_nam{loop}
fiel.(a) = fr;
fr = fr';
fr = flipud(fr);
fr(:,1:3)
f2 = fields(:,loop);
fr2 = reshape(f2,xle,yle);
fi.(field_nam{loop}) = fr2;
a = field_nam{loop};
fiel.(a) = fr2;
fr2 = fr2';
fr2 = flipud(fr2);
fr2(:,1:3)
isequal(fr,fr2)

%figure(loop)
%pcolor(x,y,fr')
%shading flat
%title(field_nam{loop})
end
%fr2     =   fr2';