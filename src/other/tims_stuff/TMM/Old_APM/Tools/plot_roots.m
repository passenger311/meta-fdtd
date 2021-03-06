function root = plot_roots(roots,contour)


if roots == 0
    FolderName = 'tmp';
    FileName    =   'Modes.txt';
    FileName    =   sprintf('%s/%s',FolderName,FileName);
    
    if exist(FileName, 'file') == 0
        FolderName  =   'Inputs';
        FileName2   =   'Modes.txt';
        FileName    =   sprintf('%s/%s',FolderName,FileName2);
        if exist(FileName, 'file') == 0
            clc
            frpintf('Modes.txt file not found')
        else
            copyfile(FileName, 'tmp');
        end
        
    end
    
    
   
    FileID      =   fopen(FileName, 'r');
    loop = 0;
    
    while 1
        loop = loop+1;
        tline = fgetl(FileID);
        if ~ischar(tline),   break,   end
        
        if loop == 6
            num = sscanf(tline(16:end),'%g');
            a = zeros(num,7);
        end
        
        if loop > 8
            a(loop-8,:) = sscanf(tline,'%g');
        end
    end
    roots = a;
    
else
    roots_tmp = zeros(size(roots,1),7);
    roots_tmp(:,2) = real(roots(:,1));
    roots_tmp(:,3) = imag(roots(:,1));
    roots_tmp(:,4) = roots(:,2);
    roots_tmp(:,5) = roots(:,3);
    roots_tmp(:,6) = roots(:,4);
    roots_tmp(:,7) = roots(:,5);
    roots = roots_tmp;
    
end

z_mat = complex(roots(:,2),roots(:,3));
ind1    =   find(roots(:,end)==1);
ind2    =   find(roots(:,end)==0);
low     =   contour(1);
high    =   contour(2);


% Plot root positions (bound roots in red, unbound in blue)
figure(1)
scatter(real(z_mat(ind1)),imag(z_mat(ind1)),30,'r')
hold on
scatter(real(z_mat(ind2)),imag(z_mat(ind2)),30,'b')
axis([low high -pi pi])
title('Roots on the complex Z plane','FontSize',12);
xlabel('Real Z','FontSize',12)
ylabel('Imag Z (pi)','FontSize',12)
set(gca,'YTick',[-pi -pi/2 0 pi/2 pi])
set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'})

if low<0
    set(gca,'XTick',[low 0 high]);
else
    set(gca,'XTick',[low high]);
end
grid on
set(gcf,'color','white');

for K = 1:length(z_mat)
    text(real(z_mat(K)),(imag(z_mat(K)))-0.2,sprintf('%.0f',K))
end
clc
view = 0;
neff_mat = roots(:,4:end);


while view~=1
    view = input('\nWould you like to view refractive index results for roots? (all/bound/unbound/n)\n','s');
    ansArr = {'n';'N';'no';'No'};
    view1 = strcmp(deblank(view),ansArr);
    view2 = strcmp(view(1),'a');
    view3 = strcmp(view(1),'b');
    view4 = strcmp(view(1),'u');
    view = max(view1);
    
    
    if view == 0
        fprintf('\n');
        for loop = 1:length(z_mat);
            
            dec     =   roots(loop,end);
            
            
            if dec == 1;
                b = sprintf('bound  ');
            else
                b = sprintf('unbound');
            end
            
            
            if sign(neff_mat(loop,2))==1
                A = sprintf('Root %g is %s with refractive index %.3g + %.3gi and loss = %.3g /cm\n'...
                    ,loop,b,neff_mat(loop,1),neff_mat(loop,2),neff_mat(loop,3));
                
            else
                A =   sprintf('Root %g is %s with refractive index %.3g - %.3gi and loss = %.3g /cm\n'...
                    ,loop,b,neff_mat(loop,1),abs(neff_mat(loop,2)),neff_mat(loop,3));
            end
            
            if(view2 == 1)
                fprintf(A);
            elseif (view3 == 1)
                if strcmp(b(1),'b')==1
                    fprintf(A);
                end
            elseif (view4 == 1)
                if strcmp(b(1),'u')==1
                    fprintf(A);
                end
            end
        end
    end
end

% Read root to track
a = sprintf('\nPlease select root(1 to %.0f/all/bound/unbound)\n',length(z_mat));
root = input(a,'s');
if strcmpi(root(1),'a')==1
    root = 1:length(z_mat);
elseif strcmpi(root(1),'b')==1
    root = ind1;
elseif strcmpi(root(1),'u')==1
    root = ind2;
else
    root = str2num(root);
end
