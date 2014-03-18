

% Sim parameters
lam = 1.55e-6;
k0 = 2*pi/lam;


% Read current job id
jfilename = 'Results/Job_name.txt';
fileID = fopen(jfilename,'r');
tline = fgetl(fileID);
name = tline;
fclose(fileID);
% Read Current velocity

vfilename = 'Results/evel_out.txt';
fileID  =   fopen(vfilename,'r');
tline   =   fgetl(fileID);

ind(1) = strfind(tline,'=');
ind(2) = strfind(tline,'/');

vel = str2double(tline(ind(1)+1:ind(2)-2));

for loop = 1:4
    tline2 = fgetl(fileID);
end

ind(1) = strfind(tline2,'=');
ind(2) = strfind(tline2,'/');
FWHM = str2double(tline2(ind(1)+1:ind(2)-2));
fclose(fileID);

% Read Current angle
a = cd;
ind1 = strfind(a,'/');
ind2 = strfind(a,'_');

angle = str2double(a(ind1(end)+1:ind2(end)-1));

% Read previous results if they exist
k = cosd(angle)*k0;
rfilename = sprintf('%sVelocity_vs_angle.txt',a(1:ind1(end)));
if exist(rfilename,'file')~=0
    in_mat = dlmread(rfilename,'\t');
    out_mat = [in_mat; angle k vel FWHM];
else
    out_mat = [angle k vel FWHM];
end
num = size(out_mat,1);
s_out_mat = sortrows(out_mat,1);
% Write results to file
fileID = fopen(rfilename,'w');
for n = 1:num
    outli = sprintf('%g\t%g\t%g\t%g\n',out_mat(n,1),out_mat(n,2),out_mat(n,3),out_mat(n,4));
    fprintf(fileID,outli);
end
fclose(fileID);

% % Plot results
% f = figure('visible','off');
% plot(s_out_mat(:,1),s_out_mat(:,3),'-o');
% xlabel('Angle (degrees)');
% ylabel('Energy velocity (c)');
% filename = sprintf('%sVelocity_angle.eps',a(1:ind1(end)));
% print(f,'-depsc','-r200','-painters',filename);
% f = figure('visible','off');
% plot(s_out_mat(:,2),s_out_mat(:,3),'-o');
% xlabel('Wavevector, \beta (m^{-1})');
% ylabel('Energy velocity (c)');
% filename = sprintf('%sVelocity_k.eps',a(1:ind1(end)));
% print(f,'-depsc','-r200','-painters',filename);
% f = figure('visible','off');
% plot(s_out_mat(:,2),s_out_mat(:,4),'-o');
% xlabel('Wavevector, \beta (m^{-1})');
% ylabel('Energy velocity (c)');
% filename = sprintf('%sFWHM_k.eps',a(1:ind1(end)));
% print(f,'-depsc','-r200','-painters',filename);



% Calculate new angle
if num == 1
    inc = -5;
    new_angle = angle+inc;
    
    % Plot results
    f = figure('visible','off');
    plot(s_out_mat(:,1),s_out_mat(:,3),'-o');
    xlabel('Angle (degrees)');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_angle.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,3),'-o');
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,4),'-o');
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sFWHM_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    
    
    
    
elseif num <= 8
    min_vel = min(s_out_mat(:,3));
    max_vel = max(s_out_mat(:,3));
    s = sign(min_vel)/sign(max_vel);
    inc = min(s_out_mat(2:end,1)-s_out_mat(1:end-1,1));
    if s == 1
        inc = -5;
        new_angle = angle+inc;
    else
        inc = inc/10;
        xnew = s_out_mat(1,1):inc:s_out_mat(end,1);
        ynew = interp1(s_out_mat(:,1),s_out_mat(:,3),xnew,'spline');
        new_angle = xnew(abs(ynew)==min(abs(ynew)));
    end
    
    
    xnew1 = s_out_mat(1,1):0.1:s_out_mat(end,1);
    ynew1 = interp1(s_out_mat(:,1),s_out_mat(:,3),xnew1,'spline');
    inc2 = (s_out_mat(end,2)-s_out_mat(1,2))/1000;
    xnew2 = s_out_mat(1,2):inc2:s_out_mat(end,2);
    ynew2 = interp1(s_out_mat(:,2),s_out_mat(:,3),xnew2,'spline');
    ynew3 = interp1(s_out_mat(:,2),s_out_mat(:,4),xnew2,'spline');
    % Plot results
    f = figure('visible','off');
    plot(s_out_mat(:,1),s_out_mat(:,3),'o',xnew1,ynew1);
    xlabel('Angle (degrees)');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_angle.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,3),'o',xnew2,ynew2);
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,4),'o',xnew2,ynew3);
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sFWHM_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    
    
    
else
    min_vel = min(s_out_mat(:,3));
    max_vel = max(s_out_mat(:,3));
    s = sign(min_vel)/sign(max_vel);
    inc = min(s_out_mat(2:end,1)-s_out_mat(1:end-1,1));
    if s == 1
        diff = sign(s_out_mat(2:end,3)-s_out_mat(1:end-1,3));
        diff2 = diff(2:end)-diff(1:end-1);
        ind = find(diff2==2);
        if isempty(ind)==1
            new_angle = max(s_out_mat(:,1))+2;
        else
            new_angle = (s_out_mat(ind(1)+2,1)+s_out_mat(ind(1)+1,1))/2;
        end
        
    else
        inc = inc/10;
        xnew = s_out_mat(1,1):inc:s_out_mat(end,1);
        ynew = interp1(s_out_mat(:,1),s_out_mat(:,3),xnew,'spline');
        new_angle = xnew(abs(ynew)==min(abs(ynew)));
    end
    
    
    xnew1 = s_out_mat(1,1):0.1:s_out_mat(end,1);
    ynew1 = interp1(s_out_mat(:,1),s_out_mat(:,3),xnew1,'spline');
    inc2 = (s_out_mat(end,2)-s_out_mat(1,2))/1000;
    xnew2 = s_out_mat(1,2):inc2:s_out_mat(end,2);
    ynew2 = interp1(s_out_mat(:,2),s_out_mat(:,3),xnew2,'spline');
    ynew3 = interp1(s_out_mat(:,2),s_out_mat(:,4),xnew2,'spline');
    % Plot results
    f = figure('visible','off');
    plot(s_out_mat(:,1),s_out_mat(:,3),'o',xnew1,ynew1);
    xlabel('Angle (degrees)');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_angle.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,3),'o',xnew2,ynew2);
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sVelocity_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    f = figure('visible','off');
    plot(s_out_mat(:,2),s_out_mat(:,4),'o',xnew2,ynew3);
    xlabel('Wavevector, \beta (m^{-1})');
    ylabel('Energy velocity (c)');
    filename = sprintf('%sFWHM_k.eps',a(1:ind1(end)));
    print(f,'-depsc','-r200','-painters',filename);
    
    
    
end





fname = sprintf('/data/users/twp11/Running/%s',name);
% Write to file
if isnan(new_angle)~=1
    thi = str2double(a(ind1(end-1)+1:ind1(end)-3));
    afilename = sprintf('%s/Angle.txt',fname);
    fileID = fopen(afilename,'w');
    ang_out = sprintf('%.3f',new_angle);
    fprintf(fileID,ang_out);
    fclose(fileID);
    mfilename = sprintf('%s/Mag.txt',fname);
    fileID = fopen(mfilename,'w');
    mag_out = sprintf('%g',thi);
    fprintf(fileID,mag_out);
    fclose(fileID);
    
else
    qfilename = sprintf('%s/quit.txt',fname);
    fileID = fopen(qfilename,'w');
    fprintf(fileID,'1');
    fclose(fileID);
end



quit;





