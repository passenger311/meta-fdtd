% Read file names
close all
clear all
clc
tstart = tic;
% Read params

filename    =   'sim_params.txt';
fileID      =   fopen(filename,'r');
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
typ         =   in{1}{1};
typ         =   isequal(typ,'dx');
if typ == 1
    resolution  =   str2double(in{1}{2});
else
    res     =   str2double(in{1}{2});
end

tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xstart      =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
inds        =   find(tline == ' ');
xfinish     =   str2double(tline(inds+1:end));
tline       =   fgetl(fileID);
in          =   textscan(tline,'%s');
central_lambda  =   str2double(in{1}{2});
fclose(fileID);


lightSI = 299792458;
%central_omega       =   1.4e15;
%central_freq        =   central_omega/(2*pi);
central_lambda      =   central_lambda*1e-6;

if  typ == 1
    resolution  =   resolution*1e-6;
else
    
    resolution          = central_lambda/res/(sqrt(11.68));
end
Courant             = 0.7;
TimeStepREAL        = (Courant*resolution)/lightSI;
time                = 0;
time                = round(time/TimeStepREAL);
k0 = 2*pi/central_lambda;


norm = 0;
%cd 28132.localhost/Data/
list = dir('Hz_*.gpl');
num = max(size(list));

for n = 1:num
    name = list(n).name;
    ind1 = strfind(name,'_');
    ind2 = strfind(name,'.0');
    str_num = str2num(name(ind1+1:ind2-1));
    num_list(n,1) = {str_num};

    num_list(n,2) = {(name(ind1+1:ind2-1))};
end

num_list = sortrows(num_list,1);
ind = [num_list{:,1}] > time;
num_list = num_list(ind,:);

num = max(size(num_list));
% cd ..
cd Results/

filename = 'fit_params.txt';
FileID      =   fopen(filename, 'r');
tline = fgetl(FileID);
loop = 0;

while 1
    loop = loop+1;
    tline = fgetl(FileID);
    if ~ischar(tline),   break, end
    
    if isempty(tline) == 0
        fit_mat(loop,:) = str2num(tline);
    end
end

amp_mat = fit_mat(:,2);
centre_mat = fit_mat(:,5);
width_mat = fit_mat(:,6);
% 
% filename = 'fit_params2.txt';
% FileID      =   fopen(filename, 'r');
% tline = fgetl(FileID);
% loop = 0;
% 
% while 1
%     loop = loop+1;
%     tline = fgetl(FileID);
%     if ~ischar(tline),   break, end
%     
%     if isempty(tline) == 0
%         fit_mat2(loop,:) = str2num(tline);
%     end
% end
% 
% amp_mat2 = fit_mat2(:,2);
% centre_mat2 = fit_mat2(:,5);
% width_mat2 = fit_mat2(:,6);
% 
% filename = 'fit_params3.txt';
% FileID      =   fopen(filename, 'r');
% tline = fgetl(FileID);
% loop = 0;
% 
% while 1
%     loop = loop+1;
%     tline = fgetl(FileID);
%     if ~ischar(tline),   break, end
%     
%     if isempty(tline) == 0
%         fit_mat3(loop,:) = str2num(tline);
%     end
% end
% 
% amp_mat3 = fit_mat3(:,2);
% centre_mat3 = fit_mat3(:,5);
% width_mat3 = fit_mat3(:,6);
% 
% filename = 'fit_params4.txt';
% FileID      =   fopen(filename, 'r');
% tline = fgetl(FileID);
% loop = 0;
% 
% while 1
%     loop = loop+1;
%     tline = fgetl(FileID);
%     if ~ischar(tline),   break, end
%     
%     if isempty(tline) == 0
%         fit_mat4(loop,:) = str2num(tline);
%     end
% end
% 
% amp_mat4 = fit_mat4(:,2);
% centre_mat4 = fit_mat4(:,5);
% width_mat4 = fit_mat4(:,6);

%cd 28132.localhost/Data/
cd ..
for n = 1:num
    clc
    fprintf('\nCycle %g out of %g',n,num)
    
    ncyc = num_list{n,2};
    filename = sprintf('Hz_%s.0.gpl',ncyc);
    file = fopen(filename,'r');
    [A,COUNTA] = textscan(file,'%n',1e5,'Delimiter','\n','CommentStyle','#');
    fclose(file);
    A = A{1};
    ind  = isnan(A) == 0;
    A = A(ind);
    
    %A = A(xstart:xfinish);
    
    A_mat(n,:) = A;
end
% cd ..
% cd ..
d   =   round(num/2);
num_list = [num_list{:,1}];
num_list = num_list';
time = num_list*TimeStepREAL;
x = 1:length(A);

for n = 1:length(amp_mat)
    a1 = amp_mat(n);
    b1 = centre_mat(n);
    c1 = width_mat(n);
%     a2 = amp_mat2(n);
%     b2 = centre_mat2(n);
%     c2 = width_mat2(n);
%     a3 = amp_mat3(n);
%     b3 = centre_mat3(n);
%     c3 = width_mat3(n);
%     a4 = amp_mat4(n);
%     b4 = centre_mat4(n);
%     c4 = width_mat4(n);    
    
    envel_mat1(n,:) = a1.*exp(-(((x-b1).^2)./(2*(c1^2))));
%     envel_mat2(n,:) = (a2.*exp(-(((x-b2).^2)./(2*(c2^2)))));
%     envel_mat3(n,:) = (a3.*exp(-(((x-b3).^2)./(2*(c3^2)))));
%     envel_mat4(n,:) = (a4.*exp(-(((x-b4).^2)./(2*(c4^2)))));
    
end

% 
% f = figure('visible','off');
% 
% 
% 
% 
% 
% 
% figure(2)
% set(gcf,'doublebuffer','on')
% sh = uicontrol('Style','slider',...
%     'Max',num,'Min',1,'Value',d,...
%     'SliderStep',[1/num 0.2],...
%     'Position',[100 0 150 30]);
% 
% 
% %the mouse-click sets a flag for point aquisition
% set(gcf,'windowbuttondownfcn','mousedown=1;');
% set(gcf,'windowbuttonupfcn','mouseup=1;');
% set(gcf,'windowbuttonmotionfcn','mousemotion=1;');
% 
% %make a control to stop the loop
% uicontrol('style','pushbutton',...
%     'string','Quit', ...
%     'position',[0 0 50 20], ...
%     'callback','stopit=1;');
% 
% %make a control to stop the loop
% hObject = uicontrol('style','checkbox',...
%     'string','Normalise', ...
%     'position',[0 20 80 20]);
% 
% %start looping and waiting for a mouse click
% stopit = 0;
% mousedown = 0;
% mouseup = 0;
% mousemotion = 0;
% ind = get(sh,'Value');
% conv    =   (resolution/central_lambda);
% x_real = x.*conv;
% 
% while (stopit==0)
%     
%     %check for valid object and chek for line 1
%     %and see if the mouse was clicked
%     indp = ind;
%     ind = get(sh,'Value');
%     
%     if indp~=ind
%         
%         %get the mouse position in graph units
%         
%         ind = get(sh,'Value');
%         ind = round(ind);
%         
%         if (get(hObject,'Value') == get(hObject,'Max'))
%             mx = 1/(max(abs(envel_mat1(ind(1),:))));
%         else
%             mx = 1/(max(abs(envel_mat1(1,:))));
%         end
%         
% 
%         plot(x_real,A_mat(ind(1),:).*mx,'LineWidth',2,'Color','b')
%         hold on
%         %plot(x_real,envel_mat1(ind(1),:).*mx,'LineWidth',2,'Color','r')
%         %hold on
% %         plot(x_real,envel_mat2(ind(1),:).*mx,'LineWidth',2,'Color','r')
% %         hold on
% %         plot(x_real,envel_mat3(ind(1),:).*mx,'LineWidth',2,'Color','r')
% %         hold on
% %         plot(x_real,envel_mat4(ind(1),:).*mx,'LineWidth',2,'Color','r')
% %         hold on
%         plot([centre_mat(ind(1)) centre_mat(ind(1))].*conv,[-2 2],'--k')
% %         hold on
% %         plot([centre_mat2(ind(1)) centre_mat2(ind(1))].*conv,[-2 2],'--k')
% %         hold on
% %         plot([centre_mat3(ind(1)) centre_mat3(ind(1))].*conv,[-2 2],'--k')
% %         hold on
% %         plot([centre_mat4(ind(1)) centre_mat4(ind(1))].*conv,[-2 2],'--k')
%         hold off
% 
%         axis([x_real(1) x_real(end) -1.1 1.1])
%         set(gcf,'color','white');
%         a   =   sprintf('Time = %5.2f fs',time(ind(1))*1e15);
%         title(a,'Fontsize',12);
%         ylabel('Hz Field, (arb)','Fontsize',12);
%         xlabel('X position (Central Wavelength)','Fontsize',12);
% 
%     end
%     
%     drawnow
% end

% Prepare the new file.
vidObj = VideoWriter('Results/peaks_norm.avi');

% Set and view the frame rate.
vidObj.FrameRate = 10;
open(vidObj);
%scrsz = get(0,'ScreenSize');
f = figure('visible','off','Position',[1 768 960 768]);

%set(gca,'nextplot','replacechildren');

for n = 2:2:num
        mx = 1;%/(max(abs(envel_mat1(n,:))));
        %mx = 1/(max(abs(envel_mat1(1,:))));
        plot(x_real,A_mat(n,:).*mx,'LineWidth',2,'Color','b')
        hold on
        plot(x_real,envel_mat1(n,:).*mx,'LineWidth',2,'Color','r')
        hold on
%         plot(x_real,envel_mat2(n,:).*mx,'LineWidth',2,'Color','r')
%         hold on
%         plot(x_real,envel_mat3(n,:).*mx,'LineWidth',2,'Color','r')
%         hold on
%         plot(x_real,envel_mat4(n,:).*mx,'LineWidth',2,'Color','r')
%         hold on
        plot([centre_mat(n) centre_mat(n)].*conv,[-2 2],'--k')
%         hold on
%         plot([centre_mat2(n) centre_mat2(n)].*conv,[-2 2],'--k')
%         hold on
%         plot([centre_mat3(n) centre_mat3(n)].*conv,[-2 2],'--k')
%         hold on
%         plot([centre_mat4(n) centre_mat4(n)].*conv,[-2 2],'--k')
        hold off

        axis([x_real(1) x_real(end) -1.1 1.1])
        set(gcf,'color','white');
        a   =   sprintf('Time = %5.2f fs',time(n)*1e15);
        title(a,'Fontsize',12);
        ylabel('Hz Field, (arb)','Fontsize',12);
        xlabel('X position (Central Wavelength)','Fontsize',12);
       set(gcf,'Renderer','painters')
           % Write each frame to the file.
   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);



% figure(3)
% ind = 480;
% int = 38; %fs
% int = int*1e-15;
% int_t = time(2)-time(1);
% int = round(int/int_t);
% x_new = (1:size(envel_mat,2)).*(resolution/central_lambda);
% a = 1/amp_mat(ind);
% plot(x_real,a.*envel_mat(ind(1),:));
% hold on
% plot(x_real,a.*envel_mat(ind(1)+int,:),'r');
% hold on
% plot(x_real,a.*envel_mat(ind(1)+2*int,:),'g');
% hold on
% plot(x_real,a.*envel_mat(ind(1)+255,:),'k');
% hold on
% plot(x_real,a.*A_mat(ind(1)+255,:),'k');

% 
% out_file = 'Results/pulses.txt';
% fileID = fopen(out_file,'w');
% fprintf(fileID,'x(lambda)\tenv_t1\tenv_t2\tenv_t3\tenv_t4\tpulse_t4\n');
% for n = 1:length(x_new)
%     x = x_real(n);
%     env_1 = a.*envel_mat(480,n);
%     env_2 = a.*envel_mat(557,n);
%     env_3 = a.*envel_mat(634,n);
%     env_4 = a.*envel_mat(735,n);
%     pulse = a.*A_mat(735,n);
%     
%     fprintf(fileID,'%g\t%g\t%g\t%g\t%g\t%g\n',x,env_1,env_2,env_3,env_4,pulse);
% end

fclose('all');

quit;