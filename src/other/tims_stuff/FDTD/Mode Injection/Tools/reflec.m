% Reflection calculation

filename{1}   =   'ebal_reflectiont.0.gpl';
filename{2}   =   'ebal_reflectioninc.0.gpl';

for loop = 1:2
    fileID  =   fopen(filename{loop},'r');
    
    for loop2   =   1:13
        tline   =   fgetl(fileID);
    end
    
    loop2   =   0;
    data1   =   [];
    
    while 1
        loop2   =   loop2+1;
        tline   =   fgetl(fileID);
        
        if ~ischar(tline),   break, end
        
        a       =   textscan(tline,'%f');
        data1(loop2,:)   =   a{1}';
               
    end
    
    energy{loop}  =   data1(:,1);
    
end
diff    =   energy{2}-energy{1};
incav   =   mean(diff(1150:1250));
refav   =   mean(energy{1}(1150:1250));
reflection  =   refav*100/incav
