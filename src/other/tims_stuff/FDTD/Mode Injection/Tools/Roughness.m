

%  close all
%  clear all
% clc
%function Roughness(noise)
noise = 1;
%Inputs
if noise == 0
    noise
    quit;
end

filename    =   'Seed.txt';
fileID      =   fopen(filename,'r');
tline = fgetl(fileID);
sd1  =   tline;
fclose(fileID);

struc_l_o     =   3000;     % Length of structure (grid cells)
dx      =   1;            % Distance between variations (grid cells)
dy      =   noise/5;       % Roughness, height variation (grid cells)
inr     =   10;             % Interpolation samples

eps1    =   1;
eps2    =   2;

ad      =   dx*4;
struc_l =   struc_l_o+ad;

list = dir('geo_*.in');
n   =   0;
for loop = 1:length(list)
    a   =   list(loop).name;
    b = find(a == 'f');
    if isempty(b) == 0
        n   =   n+1;
        ind = find(a == '_');
        nm{n}  =   a(ind(2)+1:ind(3)-1);
    end
end

nm{4} = 'prism';

for ol = 1:length(list)
    
    if ol<4
        
        filename1    =   sprintf('geo_filling_%s*',nm{ol});
        filename2   =   sprintf('%s',ls(filename1));
        sd  =   str2num(sd1(1:ol))
        
    else
        filename1   =   sprintf('geo_%s*',nm{ol-3});
        filename2   =   sprintf('%s',ls(filename1));
        sd  =   str2num(sd1(1:ol-3))
    end
    
    
    filename1 = filename2(1:end-1)
    
    %filename1    =   'geo_prism.in';
    fileID      =   fopen(filename1,'r');
    
    
    loop    =   0;
    
    while 1
        
        loop = loop+1;
        
        tline = fgetl(fileID);
        
        if isequal(tline,')SET'), break,end
        if loop == 4
            sz(1,:)    =   str2num(tline);
            info{loop} = tline;
        elseif loop > 4
            fields(loop-4,:) = str2num(tline);
        else
            info{loop} = tline;
        end
        
        
        
        
    end
    fclose(fileID);
    xle = (sz(2)-sz(1))+1
    yle = (sz(5)-sz(4))+1
    
    f = fields(:,1);
    fr = reshape(f,xle,yle);
    
    eps1    =   fr(2,1);
    eps2    =   fr(2,end);
    
    return
    struc_l =   xle+ad;
    
    %Create roughness signal
    rng(sd);
    r       =   -dy + (dy*2).*rand(struc_l/dx,1);
    size(r)
    r_int   =   interp(r,inr);
    size(r_int)
    r_int   =   r_int((inr/2)*4+2:end-(inr*4/2));
    inc     =   dx/inr;
    x0      =   0:dx:struc_l-dx-(ad);
    x       =   (inc):(inc):struc_l-(ad)-inc;
    r       =   r(4/2+1:end-4/2);
    %Find cut-throughs
    
    mx      =   max(max(abs(r_int)));
    mx      =   2*floor(mx/0.5);
    
    inc     =   1/inc;
    mx_c    =   -((mx*0.5)/2);
    r_int   =   r_int+mx_c;
    x_z_mat =   [];
    y_mat   =   [];
    
    for loop1 = 1:mx+1
        
        s_r     =   sign(r_int);
        dif     =   s_r(2:end)+s_r(1:end-1);
        ind     =   find(dif == 0);
        
        
        
        for loop = 1:length(ind)
            y1 = r_int(ind(loop));
            y2 = r_int(ind(loop)+1);
            x1 = x(ind(loop));
            x2 = x(ind(loop)+1);
            m   =   (y2-y1)/(x2-x1);
            c   =   y2-(m*x2);
            z_p =   -c/m;
            x_z(loop) = z_p;
        end
        if isempty(ind) == 1
            
        else
            x_z_mat =   [x_z_mat x_z];
            y       =   zeros(1,length(ind));
            y       =   y - mx_c;
            
            y_mat   =   [y_mat y];
        end
        
        
        r_int   =   r_int + 0.5;
        mx_c    =   mx_c + 0.5;
        x_z     =   [];
        
    end
    
    
    r_int = r_int-(((mx*0.5)/2)+0.5);
    
    x_new = [x x_z_mat];
    r_int_new = [r_int; y_mat'];
    mat = [x_new' r_int_new];
    s_mat   =   sortrows(mat,1);
    %     plot(x0,r,'k')
    %     hold on
    %     plot(s_mat(:,1),s_mat(:,2),'r')
    %     hold on
    %     scatter(x_z_mat,y_mat)
    %axis([0 10 -dy dy])
    
    r_new   =   s_mat(:,2);
    ind     =   r_new <=0;
    r_up    =   r_new;
    r_up(ind)   =   0;
    ind     =   find(r_new >= 0);
    r_lw    =   r_new;
    r_lw(ind)   =   0;
    r_lw    =   -r_lw;
    
    
    for n = 1:(mx/2)+1
        
        up{n} = r_up;
        st  =   (n-1)*0.5;
        fi  =   n*0.5;
        ind     =   up{n} <= st;
        up{n}(ind)   =   st;
        ind     =   up{n} >= fi;
        up{n}(ind)   =   fi;
        
        lw{n} = r_lw;
        st  =   (n-1)*0.5;
        fi  =   n*0.5;
        ind     =   lw{n} <= st;
        lw{n}(ind)   =   st;
        ind     =   lw{n} >= fi;
        lw{n}(ind)   =   fi;
        
        
    end
    
    
    
    
    mx_c    =   -((mx*0.5)/2);
    r_int   =   r_int+mx_c;
    x_z_mat =   [];
    y_mat   =   [];
    h_gc    =   inr*0.5/dx;
    
    for n   =   1:2
        
        if n == 1
            r_mat = up;
            epsa = eps2;
            epsb = eps1;
        else
            r_mat = lw;
            epsa = eps1;
            epsb = eps2;
        end
        
        for loop1 = 1:(mx/2)+1
            
            r   =   r_mat{loop1}-((loop1-1)*0.5);
            fi  =   0.5;
            int =   0;
            int_mat = 0;
            gc = 1;
            nx = 0;
            
            for loop2 = 1:length(s_mat)-1
                
                
                fi  =   fi + nx;
                nx  =   0;
                ind1 = loop2;
                ind2 = ind1+1;
                a   =   s_mat(ind1,1);
                b   =   s_mat(ind2,1);
                fa  =   r(ind1);
                fb  =   r(ind2);
                
                int =   int + (b-a)*((fa+fb)/2);
                
                if floor(b/fi) == 1
                    int_mat(gc) = int/(0.5^2);
                    int = 0;
                    
                    nx = 0.5;
                    av_mat(gc)  = epsb+(epsa-epsb)*int_mat(gc);
                    gc  =   gc + 1;
                    
                end
                
                
                
                
                
                
                
            end
            
            if n == 1
                av_out_up(loop1,:) = av_mat;
            else
                av_out_lw(loop1,:) = av_mat;
            end
            
        end
        
    end
    
    av_out_up = flipud(av_out_up);
    
    av_out = [av_out_up;av_out_lw];
    a   =   (mx/2)+1;
    if isequal(floor(a/2),ceil(a/2)) == 1
        n   =   1;
    else
        n   =   2;
    end
    
    
    top     =   ones(10,length(av_out))*eps1;
    bot     =   ones(10,length(av_out))*eps2;
    
    av_out = [top(1:n+8,:);av_out;bot(1:n+8,:)];
    
    av_out  =   flipud(av_out);
    %     figure(2)
    %     pcolor(av_out(2:end-2,1:20))
    %     shading interp
    av_out = [ones(size(av_out,1),1) av_out ones(size(av_out,1),3)];
    hgt     =   size(av_out,1);
    h_gc    =   floor(hgt/4);
    av_out  =   flipud(av_out);
    l = 0;
    size(av_out)
    ez_mat = [];
    ey_mat = [];
    hz_mat = [];
    hy_mat = [];
    
    for n = 1:2:hgt
        l   =   l+1;
        l2  =   0;
        
        for loop = 1:2:length(av_out)-1
            l2  =   l2+1;
            ez  =   av_out(n,loop)+av_out(n,loop+1);
            ez  =   ez + av_out(n+1,loop)+av_out(n+1,loop+1);
            ez  =   ez/4;
            
            
            ez_mat(l,l2) = ez;
            
            
        end
        l2  =   0;
        for loop = 2:2:length(av_out)-2
            l2  =   l2+1;
            hy  =   av_out(n,loop)+av_out(n,loop+1);
            hy  =   hy+av_out(n+1,loop)+av_out(n+1,loop+1);
            hy  =   hy/4;
            hy_mat(l,l2) = hy;
        end
        
        
        
    end
    %ey_mat = [ey_mat ey_mat(:,1)];
    l = 0;
    
    for n = 2:2:hgt-2
        l   =   l+1;
        l2  =   0;
        
        for loop = 1:2:length(av_out)-1
            l2  =   l2+1;
            ey  =   av_out(n,loop)+av_out(n,loop+1);
            ey  =   ey + av_out(n+1,loop)+av_out(n+1,loop+1);
            ey  =   ey/4;
            
            
            ey_mat(l,l2) = ey;
            
        end
        
        l2  =   0;
        for loop = 2:2:length(av_out)-2
            l2  =   l2+1;
            hz  =   av_out(n,loop)+av_out(n,loop+1);
            hz  =   hz+av_out(n+1,loop)+av_out(n+1,loop+1);
            hz  =   hz/4;
            hz_mat(l,l2) = hz;
        end
        
    end
    
    %hy_mat = [top(1,1:length(hy_mat)); hy_mat; bot(1,1:length(hy_mat))];
    %ez_mat = [ez_mat ez_mat(:,end)];
    %ez_mat = [top(1,1:length(ez_mat)); ez_mat; bot(1,1:length(ez_mat))];
    ez_mat  =   ez_mat(:,1:end-1);
    ey_mat  =   ey_mat(:,1:end-1);
    ey_mat  =   [ey_mat; ey_mat(end,:)];
    hz_mat  =   [hz_mat; hz_mat(end,:)];
    %hy_mat = [top(1:3,1:length(hy_mat)); hy_mat; bot(1:4,1:length(hy_mat))];
    %ez_mat = [top(1:3,1:length(hy_mat)); ez_mat; bot(1:4,1:length(hy_mat))];
    %ey_mat = [top(1:4,1:length(hy_mat)); ey_mat; bot(1:4,1:length(hy_mat))];
    %hz_mat = [top(1:4,1:length(hy_mat)); hz_mat; bot(1:4,1:length(hy_mat))];
    %     figure(3)
    %     pcolor(hz_mat)
    %     shading interp
    %     figure(4)
    %     pcolor(ey_mat)
    %     shading interp
    %     figure(5)
    %     pcolor(hy_mat)
    %     shading interp
    %     figure(6)
    %     pcolor(ez_mat)
    %     shading interp
    
    ez_mat  =   ez_mat';
    ey_mat  =   ey_mat';
    hz_mat  =   hz_mat';
    hy_mat  =   hy_mat';
    
    out.ex  =   reshape(hy_mat,xle*yle,1);
    out.ey  =   reshape(ey_mat,xle*yle,1);
    out.ez  =   reshape(ez_mat,xle*yle,1);
    out.hx  =   reshape(ey_mat,xle*yle,1);
    out.hy  =   reshape(hy_mat,xle*yle,1);
    out.hz  =   reshape(hz_mat,xle*yle,1);
    
    filename    =   filename1;
    fileID      =   fopen(filename,'w');
    
    for loop = 1:length(info)
        fprintf(fileID,'%s\n',info{loop});
    end
    
    for loop = 1:length(out.ex)
        fprintf(fileID,'  %g %g %g %g %g %g\n',out.ex(loop),out.ey(loop),out.ez(loop),out.hx(loop),out.hy(loop),out.hz(loop));
    end
    fprintf(fileID,')SET');
    fclose(fileID);
    close all
end


%end
