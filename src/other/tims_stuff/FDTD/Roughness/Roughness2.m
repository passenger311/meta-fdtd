%ROUGHNESS Create surface roughness
%   SURF = Roughness()....
clear all
close all
%function Roughness2(xst,xfi,yst,x_noise,y_noise,noise_type,dx,res_type,lambda,n_max)

% Read parameters in
filename    =   'R_params.txt';
fileID  =   fopen(filename,'r');
loop    =   0;
while 1
    tline   =   fgetl(fileID);
    loop    =   loop+1;
    if ischar(tline) == 0
        break
    end
    inf(loop)     =   textscan(tline,'%s');
end
fclose(fileID);
lambda  =   str2double(inf{8}{end});
xst     =   0;
xfi     =   str2double(inf{1}{end});
x_noise_min     =   str2double(inf{2}{end});
x_noise_max     =   str2double(inf{3}{end});
y_noise     =   str2double(inf{4}{end});


x_type  =   str2double(inf{5}{end});
y_type  =   str2double(inf{6}{end});
n_max   =   str2double(inf{7}{end});
res_type    =   inf{9}{1};
res     =   str2double(inf{9}{end});
x_width     =   str2double(inf{10}{end});
real.lambda         =   lambda*1e-6;
comp.xst            =   xst;
comp.xfi            =   xfi;

real.noise.xmin     =   x_noise_min*1e-9;
real.noise.xmax     =   x_noise_max*1e-9;
real.noise.y        =   y_noise*1e-9;
real.noise.x        =   x_width*1e-9;

% Calculate dx
if res_type(1) == 'r'
    real.dx     =   real.lambda/n_max/res;
else
    real.dx     =   res*1e-6;
end

real.xst            =   comp.xst*real.dx;
real.xfi            =   comp.xfi*real.dx;



% Read in interfaces

filename    =   'R_interfaces.txt';
fileID  =   fopen(filename,'r');
loop    =   0;
inf     =   [];
while 1
    tline   =   fgetl(fileID);
    loop    =   loop+1;
    if ischar(tline) == 0
        break
    end
    inf_int{loop}     =   textscan(tline,'%s\t%s\t%s\t%s\t%s');
end
fclose(fileID);

% Read seed

filename    =   'Seed.txt';
seedlist    =   dlmread(filename);
% filename2 = 'Num.txt';
% runnum = dlmread(filename2);
runnum = 1;
%seed = num2str(seedlist(runnum));
%seed    =   '6548952';
%fclose(fileID);
seed = seedlist;

num_int     =   length(inf_int);



for ml = 1:num_int
    int_name    =   inf_int{ml}{1}{1};
    yst     =   str2double(inf_int{ml}{2}{1});
    layer_num   =   str2double(inf_int{ml}{3}{1});
    epsa    =   str2num(inf_int{ml}{4}{1});
    epsb    =   str2num(inf_int{ml}{5}{1});
    comp.yst            =   yst;
    real.yst            =   comp.yst*real.dx;
    sd  =   seed(layer_num)  % str2double(seed(1:layer_num))
    rng(sd)
    real.x  =   [];
    comp.x  =   [];
    comp.r  =   [];
    r   =   [];
    x   =   [];
    comp.y  =   [];
    
    % Create spacing along x direction
    
    if x_type == 1
        real.x  =   real.xst:real.noise.x:real.xfi;
        if real.x(end) ~= real.xfi
            real.x(end+1)   =   real.xfi;
        end
    elseif x_type == 2
        a   =   real.noise.xmin;
        b   =   real.noise.xmax;
        loop    =   0;
        x   =   0;
        real.x  =   0;
        while x<real.xfi
            r   =   a + (b-a).*rand;
            x   =   x + r;
            real.x  =   [real.x x];
        end
        
        real.x(end)     =   real.xfi;
    elseif x_type == 3
        a   =   (real.noise.xmax+real.noise.xmin)/2;
        b   =   (real.noise.xmax-real.noise.xmin)/2;
        b   =   b*0.25;
        x   =   0;
        real.x  =   0;
        while x<real.xfi
            r   =   a + (b-a).*rand;
            x   =   x + r;
            real.x  =   [real.x x];
        end
        
    end
    
    l   =   length(real.x);
    % Create surface roughness
    if y_type == 1
        r   =   zeros(1,l);
        r   =   r';
    elseif y_type == 2                  % Uniformly distributed noise
        a   =   real.noise.y;
        b   =   -a;
        r   =   a + (b-a).*rand(l,1);
    elseif y_type == 3                                % Gaussian distributed noise
        a   =   real.noise.y/2;
        r   =   a.*randn(l,1);
        ind     =   r>real.noise.y;
        r(ind)  =   real.noise.y;
        ind     =   r<-real.noise.y;
        r(ind)  =   -real.noise.y;
    elseif y_type == 4
        fnametmp = 'wavenumber.txt';
        k = dlmread(fnametmp,'\t');
	fprintf('k = %g\n',k);
        a = real.noise.y;
        r = a*sin(k*(real.x + sd*100));
        r = r';
 	fprintf('max = %g\n',max(r));
    end

    r(1:60) = zeros(60,1);
    r(end-59:end) = zeros(60,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot(real.x,r)
    %r = a*sin(2*pi*(2/1.55e-6)*real.x+rand(1));
    %r = r'.*rand(l,1);
%     plot(real.x*1e6,r*1e9,'LineWidth',3)
%     set(gcf,'color','w')
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     xlabel('X position, x (\mum)',...
%         'FontSize',20,'FontName','Arial')
%     yl = sprintf('Surface Roughness, y, (nm)');
%     ylabel(yl,'FontSize',20,'FontName','Arial')
%     set(gca,'fontsize',20,'fontname','arial')
%     xlim([0 5])
%     zeroPadFactor = nextpow2(length(r)) + 3;
%     
%     Fs      =   1/(real.x(2)-real.x(1));
%     L       =   length(r);
%     NFFT    =   2^nextpow2(L);
%     Y       =   fft(r,NFFT)/L;
%     f       =   Fs/2*linspace(0,1,NFFT/2+1);
%     
%     N       =   2^zeroPadFactor;
%     k       =   0:N-1;                          %create a vector from 0 to N-1
%     T       =   N/Fs;                           %get the frequency interval
%     freq    =   k/T;                            %create the frequency range
%     X       =   fft(r,N)/length(real.x);      % normalize the data
%     
%     %only want the first half of the FFT, since it is redundant
%     cutOff  =   ceil(N/2);
%     
%     %take only the first half of the spectrum
%     X       =   X(1:cutOff);
%     freq    =   freq(1:cutOff);
%     f2      =   Fs/2*linspace(0,1,N/2+1);
%     f2      =   1./f2(1:end-1);
%     amp     =   2*abs(X(1:N/2));
%     [pks,locs] = findpeaks(amp,'MINPEAKHEIGHT',max(amp)/10,'MINPEAKDISTANCE',100);
%     peak = f2(locs)
%     figure(2)
%     semilogx(f2,amp,'LineWidth',3)
%     set(gcf,'color','w')
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     xlabel('Wavelength, (m)',...
%         'FontSize',20,'FontName','Arial')
%     yl = sprintf('Surface Roughness Spectrum, y');
%     ylabel(yl,'FontSize',20,'FontName','Arial')
%     set(gca,'fontsize',20,'fontname','arial')
%     xlim([min(f2) max(f2)])
%     text(peak*2, pks*2/3,sprintf('%g',peak),'FontSize',20,'FontName','Arial')
%     RMS = sqrt(sum(r.^2)/length(r))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    % Convert to computational units
    comp.x  =   real.x/real.dx;
    comp.r  =   r'/real.dx;
    
    %figure(2)
    %plot(comp.x,comp.r)
    %hold on
    
    % Calculate values at x grid boundarys
    x   =   0;
    x_mat   =   [];
    y_mat   =   [];
    loop    =   0;
    while x<comp.xfi
        a   =   find(comp.x==x, 1);
        if isempty(a) == 1
            loop    =   loop+1;
            ind1    =   find(comp.x<x);
            x_min   =   comp.x(ind1(end));
            ind2    =   find(comp.x>x);
            x_max   =   comp.x(ind2(1));
            y_min   =   comp.r(ind1(end));
            y_max   =   comp.r(ind2(1));
            m       =   (y_max-y_min)/(x_max-x_min);
            c       =   y_min-m*x_min;
            y_new   =   (x*m)+c;
            y_mat(loop)     =   y_new;
            x_mat(loop)     =   x;
        end
        
        x   =   x+0.5;
    end
    
    comp.x  =   [comp.x x_mat];
    comp.y  =   [comp.r y_mat];
    
    comp.mat    =   [comp.x' comp.y'];
    
    comp.mat    =   sortrows(comp.mat,1);
    comp.x      =   comp.mat(:,1);
    comp.y      =   comp.mat(:,2);
    
    comp.y  =   comp.y+comp.yst;
    
    a   =   sprintf('Layer name is %s\n',int_name);
    fprintf(1,a);
    a   =   sprintf('Interface position = %g\n',comp.yst);
    fprintf(1,a);
    a   =   sprintf('Noise position = min(%g) mean(%g) max(%g) rms(%g)\n\n',min(comp.y),mean(comp.y),max(comp.y),sqrt(mean((comp.y-mean(comp.y)).^2)));
    fprintf(1,a);
   
    % Calculate values at y grid boundarys
    
    min_y   =   min(comp.y);
    max_y   =   max(comp.y);
    mul     =   floor(abs(min_y)/0.5);
    yst     =   sign(min_y)*mul*0.5;
    mul     =   ceil(abs(max_y/0.5));
    yfi     =   sign(max_y)*mul*0.5;
    loop1   =   0;
    x_mat   =   [];
    y_mat   =   [];
    y_tmp   =   [yst yfi];
    yst     =   min(y_tmp);
    yfi     =   max(y_tmp);
    
    for y = yst:0.5:yfi
        
        tmp.y   =   sign(comp.y-y);
        diff    =   tmp.y(2:end)-tmp.y(1:end-1);
        ind     =   find(abs(diff) == 2);
        
        for loop = 1:length(ind)
            a   =   ind(loop);
            min_x   =   comp.x(a);
            max_x   =   comp.x(a+1);
            min_y   =   comp.y(a);
            max_y   =   comp.y(a+1);
            m   =   (max_y-min_y)/(max_x-min_x);
            c   =   (max_y-m*max_x);
            x_new   =   (y-c)/m;
            x_mat(loop1+loop)    =   x_new;
            y_mat(loop1+loop)    =   y;
        end
        loop1   =   length(x_mat);
        
        
    end
    
    
    comp.y  =   [comp.y;y_mat'];
    comp.x  =   [comp.x;x_mat'];
    comp.mat    =   [comp.x comp.y];
    
    comp.mat    =   sortrows(comp.mat,1);
    comp.x  =   comp.mat(:,1);
    comp.y  =   comp.mat(:,2);
    %figure(3)
    %plot(comp.x,comp.y,'r')
    xst     =   0;
    xfi     =   comp.x(end);
    y_ind   =   0;
    grid_out    =   [];
    % Calculate grid
    min_y   =   min(comp.y);
    mul     =   floor(min_y/0.5);
    yst     =   mul*0.5;
    max_y   =   max(comp.y);
    mul     =   floor(max_y/0.5);
    yfi     =   mul*0.5;
    x_mat   =   [];
    for y = yst:0.5:yfi
        y_ind   =   y_ind+1;
        tmp.y   =   comp.y-y;
        ind     =   tmp.y>0.5;
        tmp.y(ind)  =   0.5;
        ind     =   tmp.y<0;
        tmp.y(ind)  =   0;
        ind     =   1;
        x_ind   =   0;
        for x = xst:0.5:xfi-0.5
            x_ind   =   x_ind+1;
            t_int   =   0;
            while 1
                
                x_max   =   comp.x(ind+1);
                x_min   =   comp.x(ind);
                y_max   =   tmp.y(ind+1);
                y_min   =   tmp.y(ind);
                int     =   (x_max-x_min)*((y_min+y_max)/2);
                t_int   =   t_int+int;
                ind     =   ind+1;
                x_max   =   round(x_max*1e8);
                x_max   =   x_max/1e8;
                if x_max == x+0.5
                    break
                end
                
            end
           
            x_mat(x_ind) = t_int/(0.25);
        end
        grid_out(y_ind,:)   =   x_mat;
    end
    hgt     =   size(grid_out,1);
    b_pad   =   1;
    a   =   sprintf('Minimum grid value = %g\n',min(grid_out(1,:)));
    fprintf(1,a);    
    a   =   sprintf('Maximum grid value = %g\n\n',max(grid_out(y_ind,:)));
    fprintf(1,a);
    

    
    
    if round(yst) == yst;
        t_pad   =   0;
        b_pad   =   0;
    else
        t_pad   =   1;
        b_pad   =   1;
    end
    
    if yst  == yfi
        if round(yst) == yst;
            t_pad   =   1;
            b_pad   =   0;
        else
            t_pad   =   0;
            b_pad   =   1;
        end
        if yst == 0
            grid_out    =   [];
        end
    end
    
    
    filename    =   sprintf('geo_%s.in',int_name);
    fileID  =   fopen(filename,'r');
    
    
    loop    =   0;
    
    for loop = 1:4
        tline = fgetl(fileID);
        if loop == 4
            sz(1,:)    =   str2num(tline);
            inf{loop} = tline;
        else
            inf{loop} = tline;
        end
    end
    
    fclose(fileID);
    xle = (sz(2)-sz(1))+1;
    yle = (sz(5)-sz(4))+1;
    t_h     =   (sz(5)-(yfi+0.5)+1)*2;
    b_h     =   (yst-sz(4)+1)*2;
    
    
    if isempty(grid_out) == 0
        grid_out    =   flipud(grid_out);
        l   =   length(grid_out);
        s   =   1;
    else
        l   =   (2*comp.xfi);
        t_h     =   t_h+1;
    end
    
    
    
    top     =   zeros(t_h,l);
    bot     =   ones(b_h,l);
    
    grid_out    =   [top;grid_out;bot];
    if int_name(1) == 'f'
        left    =   zeros(size(grid_out,1),1);
    else
        
        left    =   ones(size(grid_out,1),1);
    end
    
    l   =   length(grid_out);
    hgt     =   size(grid_out,1);
    
    % Calculate fields
    
    grid_out    =   epsb+(epsa-epsb)*grid_out;
   
    grid_out    =   [left grid_out left];
    ez_mat = [];
    ey_mat = [];
    hz_mat = [];
    hy_mat = [];
    l   =   0;
    
    for n = 2:2:hgt-1
        l   =   l+1;
        l2  =   0;
        
        for loop = 1:2:length(grid_out)-1
            l2  =   l2+1;
            ez  =   grid_out(n,loop)+grid_out(n,loop+1);
            ez  =   ez + grid_out(n+1,loop)+grid_out(n+1,loop+1);
            ez  =   ez/4;
            
            
            ez_mat(l,l2) = ez;
            
            
        end
        l2  =   0;
        for loop = 2:2:length(grid_out)-2
            l2  =   l2+1;
            hy  =   grid_out(n,loop)+grid_out(n,loop+1);
            hy  =   hy+grid_out(n+1,loop)+grid_out(n+1,loop+1);
            hy  =   hy/4;
            hy_mat(l,l2) = hy;
        end
        
        
        
    end
    
    l = 0;
    
    for n = 1:2:hgt-2
        l   =   l+1;
        l2  =   0;
        
        for loop = 1:2:length(grid_out)-1
            l2  =   l2+1;
            ey  =   grid_out(n,loop)+grid_out(n,loop+1);
            ey  =   ey + grid_out(n+1,loop)+grid_out(n+1,loop+1);
            ey  =   ey/4;
            
            
            ey_mat(l,l2) = ey;
            
        end
        
        l2  =   0;
        for loop = 2:2:length(grid_out)-2
            l2  =   l2+1;
            hz  =   grid_out(n,loop)+grid_out(n,loop+1);
            hz  =   hz+grid_out(n+1,loop)+grid_out(n+1,loop+1);
            hz  =   hz/4;
            hz_mat(l,l2) = hz;
        end
        
    end
    
    if int_name(1) == 'f'
        hz_mat  =   [hz_mat zeros(size(hz_mat,1),1)];
        hy_mat  =   [hy_mat zeros(size(hy_mat,1),1)];
    else
        hz_mat  =   [hz_mat ones(size(hz_mat,1),1)];
        hy_mat  =   [hy_mat ones(size(hy_mat,1),1)];
    end

    
    hz  =   flipud(hz_mat);
    ez  =   flipud(ez_mat);
    ey  =   flipud(ey_mat);
    hy  =   flipud(hy_mat);
    
    filename    =    sprintf('geo_%s.in',int_name);
    fileID  =   fopen(filename,'w');
    for loop = 1:4
        fprintf(fileID,'%s\n',inf{loop});
    end
    hz  =   hz';
    hz  =   reshape(hz,xle*yle,1);
    hy  =   hy';
    hy  =   reshape(hy,xle*yle,1);
    ey  =   ey';
    ey  =   reshape(ey,xle*yle,1);
    ez  =   ez';
    ez  =   reshape(ez,xle*yle,1);
    for loop = 1:length(hz)
        fprintf(fileID,'  %g %g %g %g %g %g\n',hy(loop),ey(loop),ez(loop),ey(loop),hy(loop),hz(loop));
    end
    fprintf(fileID,')SET');
    fclose(fileID);
end
quit;
