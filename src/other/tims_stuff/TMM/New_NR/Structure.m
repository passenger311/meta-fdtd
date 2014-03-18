classdef Structure < handle
    
    properties (Constant)
        siC = 299792458;            % Speed of light (m/s)
        eps0 = 8.85418782e-12;      % Permittivity of free space
        mu0 = 1.25663706e-6;        % Permeabillity of free space
    end %Constant properties
    
    properties
        layers =[];                 % Layer definitions
        d;                          % Thickness of layers (m)
        omega = 1;                  % Free space angular frequency
        pola;                       % Mode Polarization
        nrs = 100;                  % Number of Newton Raphson steps
        convcheck = 1;              % Check Newton Raphson converges
        bnd = 0;
        beta = 1e6;
    end %Structure properties
    
    properties (SetAccess = private)
        eps;                        % Permittivity of layers
        deps;                       % Derivative of eps wrt omega
        mu;                         % Permeability of layers
        dmu;                        % Derivative of mu wrt omega
        n;                          % Refractive index of layers
        dn;                         % derivative of n wrt omega
        impz;                       % Impedance of layers
        num;                        % Length
        width;                      % Total width
        roots_w;                    % Complex omega roots
        roots_k;                    % Complex beta roots
    end %Private properties
    
    methods
        function obj = Structure(layers)
            obj.layers = layers;
            if ~isempty(obj.layers(1))
                leng = length(obj.layers);
                fid = fopen('temp.m','wt');
                fprintf(fid,'function output = temp(w)\n');
                fprintf(fid,'l = %g;\n', leng);
                
                for loop = 1:leng
                    
                    if strcmpi(obj.layers{loop}(1:2),'me')
                        obj.layers{loop}(7:end-1);
                        params = str2num(obj.layers{loop}(7:end-1));
                        tmp = obj.Metal(params(1),params(2),params(3),obj.omega,fid,loop);
                    elseif strcmpi(obj.layers{loop}(1:2),'di')
                        params = str2num(obj.layers{loop}(12:end-1));
                        tmp = obj.Dielectric(params(1),params(2),fid,loop);
                    elseif strcmpi(obj.layers{loop}(1:2),'nr')
                        params = str2num(obj.layers{loop}(5:end-1));
                        tmp = obj.NRI(params(1),params(2),params(3),obj.omega,fid,loop);
                    elseif strcmpi(obj.layers{loop}(1:2),'lo')
                        params = str2num(obj.layers{loop}(9:end-1));
                        tmp = obj.Lorentz(params(1),params(2),params(3),params(4),obj.omega,fid,loop);
                    elseif strcmpi(obj.layers{loop}(1:2),'fo')
                        params = str2num(obj.layers{loop}(10:end-1));
                        tmp = obj.four_lvl(obj,params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8),params(9),fid,loop);
                    else
                        layername = sprintf('obj.%s(%g,%g,%g)',obj.layers{loop},obj.omega,fid,loop);
                        p = eval(layername);
                        tmp = p(obj.omega);
                    end
                end
                fprintf(fid,'output = {eps,mu,deps,dmu};\n');
                fprintf(fid,'end\n');
                
            end
        end
        
        function value = get.num(obj)
            leng = length(obj.layers);
            tmp = temp(obj.omega);
            obj.eps = tmp{1};
            obj.mu = tmp{2};
            obj.deps = tmp{3};
            obj.dmu = tmp{4};
            
            
            obj.n = sqrt(obj.eps.*obj.mu);
            obj.dn = (obj.deps.*obj.mu + ...
                obj.eps.*obj.dmu)./(2*obj.n);
            obj.impz = sqrt(obj.mu./obj.eps);
            obj.width = sum(obj.d);
            
            value = leng;
        end
        
        function roots = find_roots_w(obj,wst,wfi,wres,b_in)
            w_inc = (wfi-wst)/wres;
            w_in = wst:w_inc:wfi;
            
            old = obj.nrs;
            obj.nrs = old*2;
            wg_mat = zeros(length(w_in),2);
            for loop = 1:length(w_in)
                w   =   w_in(loop);
                obj.omega = w;
                obj.num;
                w   =   obj.NR(w,b_in,1,'obj.eigenfunctionW');
                wg_mat(loop,:)    =    w;
            end
            obj.nrs = old;
            
            rts = wg_mat(~isnan(wg_mat(:,1)));
            idp = 7;
            rrts = obj.roundd(rts,idp);
            unrts = unique(rrts);
            
            for lp = 1:length(unrts)
                ind = find(rrts == unrts(lp));
                rtsnew(lp) = rts(ind(1));
            end
            
            rts = rtsnew(real(rtsnew)>wst & real(rtsnew)<wfi);
            
            
            for i = 1:length(rts)
                tmp = obj.field_plotter(rts(i),b_in,0e-9,0e-9,1000,0);
                fz = real(tmp{1}(:,2));
                ynew = sign(fz);
                ynext = ynew(2:end)+ynew(1:end-1);
                n(i) = length(find(ynext == 1 | ynext == 0));
            end
            
            roots = sortrows([n.' rts.']);
            obj.roots_w = roots(:,2).';
            
        end
        
        function roots = find_roots_k(obj,bst,bfi,bres,w_in)
            b_inc = (bfi-bst)/bres;
            b_in = bst:b_inc:bfi;
            obj.omega = w_in;
            obj.num;
            old = obj.nrs;
            obj.nrs = old*2;
            b_mat = zeros(length(b_in),2);
            
            for loop = 1:length(b_in)
                b   =   b_in(loop);
                b   =   obj.NR(w_in,b,0,'obj.eigenfunctionK');
                b_mat(loop,:)    =    b;
            end
            
            obj.nrs = old;
            
            
            rts = b_mat(~isnan(b_mat(:,1)));
            idp = 4;
            rrts = obj.roundd(rts,idp);
            unrts = unique(rrts);
            
            for lp = 1:length(unrts)
                ind = find(rrts == unrts(lp));
                rtsnew(lp) = rts(ind(1));
            end
            
            rts = rtsnew(real(rtsnew)>bst & real(rtsnew)<bfi);
            
            
            for i = 1:length(rts)
                tmp = obj.field_plotter(w_in,rts(i),0e-9,0e-9,1000,0);
                fz = real(tmp{1}(:,2));
                ynew = sign(fz);
                ynext = ynew(2:end)+ynew(1:end-1);
                n(i) = length(find(ynext == 1 | ynext == 0));
            end
            
            %             tmp = b_mat(b_mat(:,2)==1,1);
            %             tmp(isnan(real(tmp)))=[];
            %             b_mat_r = (obj.roundd(tmp,5));
            %             b_mat_r(isnan(imag(b_mat_r)))=real(b_mat_r(isnan(imag(b_mat_r))));
            %             b_mat_ru = unique(b_mat_r);
            %             b_mat = [];
            %             for loop = 1:length(b_mat_ru)
            %                 ind = find(b_mat_r==b_mat_ru(loop));
            %                 b_mat(loop) = tmp(ind(1));
            %             end
            %             b_mat = b_mat.*sign(real(b_mat));%b_mat(sign(real(b_mat))~=-1);
            %             %b_mat = b_mat(sign(imag(b_mat))==1);
            %             b_mat = b_mat(real(b_mat)<bfi);
            %             b_mat = sortrows(b_mat);
            
            roots = sortrows([n.' rts.']);
            obj.roots_k = roots(:,2).';
            
        end
        
        function mode = dispersionW(obj,wst,b_in,bst,bfi,bres)
            b_inc = (bfi-bst)/bres;
            b_max = bfi;
            data = [];
            nind2 = 0;
            perold = 0;
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            for loop = 1:2
                w = wst;
                nind = 0;
                tmp = b_in:b_inc:b_max;
                for b = b_in:b_inc:b_max
                    obj.omega = w;
                    obj.num = obj.num;
                    nind = nind + 1;
                    nind2 = nind2 + 1;
                    per = floor(1000*nind2/bres);
                    tmp = obj.NR(w,b,1,'obj.eigenfunctionW');
                    w = tmp(1);
                    gvel = obj.vg(w,b);
                    %evel = obj.ve(w,b,0,1000e-9,10000);
                    n2_val      =   obj.n.^2;
                    gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
                    gammaN      =   sgn(2)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
                    data_half(nind,:) = [b tmp(2) real(w) imag(w) gvel];
                end
                data = [data;data_half];
                b_max = bst;
                b_inc = -b_inc;
                b_in = b_in+b_inc;
                data_half = [];
            end
            
            mode = sortrows(data,1);
            
        end
        
        function mode = geodispersionW(obj,win,b,dst,dfi,dres,lay)
            d_inc = (dfi-dst)/dres;
            d_max = dfi;
            data = [];
            nind2 = 0;
            perold = 0;
            d_in = dst;
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            for loop = 1:2
                w = win;
                nind = 0;
                %tmp = d_in:d_inc:d_max;
                for thick = d_in:d_inc:d_max
                    nind = nind + 1;
                    nind2 = nind2 + 1;
                    %per = floor(1000*nind2/bres);
                    %                     if per>perold
                    %                         fprintf(1,sprintf('Calculating %g complete\n',per/10));
                    %                         perold = per;
                    %                     end
                    obj.d(lay) = thick;
                    obj.num = obj.num;
                    tmp = obj.NR(w,b,1,'obj.eigenfunctionW');
                    w = tmp(1);
                    gvel = obj.vg(w,b);
                    n2_val      =   obj.n.^2;
                    gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
                    gammaN      =   sgn(2)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
                    data_half(nind,:) = [thick tmp(2) real(w) imag(w) 2*imag(w) gvel,...
                        real(gamma1) imag(gamma1) real(gammaN) imag(gammaN)];
                    
                    %                     if nind > 2
                    %                         rew = spline(data_half(:,1),data_half(:,3),thick+d_inc);
                    %                         imw = spline(data_half(:,1),data_half(:,4),thick+d_inc);
                    %                         w = complex(rew,imw);
                    %                     end
                    
                    %
                    %                     if b <= real(w)/obj.siC
                    %                         break
                    %                     end
                end
                data = [data;data_half];
                d_max = dst;
                d_inc = -d_inc;
                d_in = d_in+d_inc;
                data_half = [];
            end
            
            mode = sortrows(data,1);
            
        end
        
        function mode = dispersionK(obj,bst,w_in,wst,wfi,wres)
            w_inc = (wfi-wst)/wres;
            w_max = wfi;
            data = [];
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            for loop = 1:2
                b = bst;
                nind = 0;
                for w = w_in:w_inc:w_max
                    obj.omega = w;
                    obj.num = obj.num;
                    nind = nind + 1;
                    tmp = obj.NR(w,b,1,'obj.eigenfunctionK');
                    b = tmp(1);
                    gvel = obj.vgk(w,b);
                    %evel = obj.ve(w,b,1000e-9,1000e-9,10000);
                    %n2_val      =   obj.n.^2;
                    %gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
                    %gammaN      =   sgn(2)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
                    data_half(nind,:) = [w tmp(2) real(b) imag(b) gvel];
                end
                data = [data;data_half];
                w_max = wst;
                w_inc = -w_inc;
                data_half = [];
            end
            
            mode = sortrows(data,1);
            
        end
        
        function mode = realdispersion(obj,bst,bfi,bres,wst,wfi,wres)
            winc = (wfi-wst)/wres;
            binc = (bfi-bst)/bres;
            mat = zeros(wres,bres);
            y = 0;
            for w = wst:winc:wfi
                y = y + 1;
                x = 0;
                obj.omega = complex(w,0);
                obj.num = obj.num;
                for b = bst:binc:bfi
                    disp = obj.eigenfunctionW(obj,complex(w,0),complex(b,0));
                    x = x + 1;
                    mat(y,x) = disp(1);
                end
                
            end
            
            mode = {bst:binc:bfi,wst:winc:wfi,mat};
            
            
        end
        
        function out = findvp(obj,b,w,vel)
            options=optimset('TolFun',1e-30,'display','off');       % Option to display output
            [x,fval] = fsolve(@opt,b,options);  % Call optimizer
            b0 = x;
            tmp = obj.NR(w,b0,1,'obj.eigenfunctionW');
            w = tmp(1);
            new = real(obj.vg(w,b0));
            out = [b0 w new];
            
            function value = opt(b)
                tmp = obj.NR(w,b,1,'obj.eigenfunctionW');
                w = tmp(1);
                value = real(obj.vg(w,b)-vel);
            end
        end
        
        function out = findb(obj,b,w,w0)
            options=optimset('TolFun',1e-30,'display','off');       % Option to display output
            [x,fval] = fsolve(@opt,b,options);  % Call optimizer
            b0 = x;
            tmp = obj.NR(w,b0,1,'obj.eigenfunctionW');
            w = tmp(1);
            new = real(obj.vg(w,b0));
            out = [b0 w new abs(w-w0)];
            
            function value = opt(b)
                tmp = obj.NR(w,b,1,'obj.eigenfunctionW');
                w = tmp(1);
                value = real(real(w)-w0);
            end
        end
        
        function out = plotDispersionRelation(obj,w,b)
            obj.omega = w;
            obj.num = obj.num;
            out = obj.eigenfunctionK(obj,w,b);
        end
        
        function out = reflec1D(obj,st,fi,res,var)
            inc = (fi-st)/res;
            loop = 0;
            switch var
                case 'Omega'
                    for w = st:inc:fi
                        loop = loop + 1;
                        obj.omega = w;
                        obj.num = obj.num;
                        val(loop,:) = obj.RTA(obj,w,obj.beta);
                    end
                    
                case 'Beta'
                    for b = st:inc:fi
                        loop = loop + 1;
                        val(loop,:) = obj.RTA(obj,obj.omega,b);
                    end
                case 'Angle'
                    for theta = st:inc:fi
                        loop = loop + 1;
                        b = (obj.omega/obj.siC)*sind(theta);
                        val(loop,:) = obj.RTA(obj,obj.omega,b);
                    end
                otherwise
                    warning('Unexpected plot type.');
                    
            end
            out = [(st:inc:fi).' val];
        end
        
        function out = reflec2D(obj,wst,wfi,wres,st,fi,res,var)
            winc = (wfi-wst)/wres;
            inc = (fi-st)/res;
            y = 0;
            TM_R = zeros(wres,res);
            TM_T = TM_R;
            TM_A = TM_R;
            TE_R = TM_R;
            TE_T = TM_R;
            TE_A = TM_R;
            
            for w = wst:winc:wfi
                obj.omega = w;
                obj.num = obj.num;
                x = 0;
                y = y + 1;
                switch var
                    
                    case 'Beta'
                        for b = st:inc:fi
                            x = x + 1;
                            tmp = obj.RTA(obj,w,b);
                            TM_R(y,x) = tmp(1);
                            TM_T(y,x) = tmp(2);
                            TM_A(y,x) = tmp(3);
                            TE_R(y,x) = tmp(4);
                            TE_T(y,x) = tmp(5);
                            TE_A(y,x) = tmp(6);
                        end
                    case 'Angle'
                        for theta = st:inc:fi
                            x = x + 1;
                            b = (obj.omega/obj.siC)*sind(theta);
                            tmp = obj.RTA(obj,w,b);
                            TM_R(y,x) = tmp(1);
                            TM_T(y,x) = tmp(2);
                            TM_A(y,x) = tmp(3);
                            TE_R(y,x) = tmp(4);
                            TE_T(y,x) = tmp(5);
                            TE_A(y,x) = tmp(6);
                        end
                    otherwise
                        warning('Unexpected plot type.');
                        
                end
                
            end

            out = {st:inc:fi,wst:winc:wfi,TM_R,TM_T,TM_A,...
                TE_R,TE_T,TE_A};
        end
        
        function out = NR(obj,w,b,pl,compmode)
            
            
            z_mat = zeros(1,obj.nrs);
            
            if strcmpi(compmode(end),'W')
                z = w;
                % Perform Newton Raphson method for nrs number of times
                for loop = 1:obj.nrs
                    obj.omega = z;
                    obj.num = obj.num;
                    sol = obj.eigenfunctionW(obj,z,b);
                    H = sol(1);
                    DH = sol(2);
                    z = z - (H/DH);     % Newton Raphson equation
                    z_mat(loop) = z;
                end
            else
                z = b;
                % Perform Newton Raphson method for nrs number of times
                for loop = 1:obj.nrs
                    sol = obj.eigenfunctionK(obj,w,z);
                    H = sol(1);
                    DH = sol(2);
                    z = z - (H/DH);     % Newton Raphson equation
                    z_mat(loop) = z;
                end
            end
            
            
            
            % Check that root converges, return 0 if not
            if pl == 1
                diff = z_mat(2:end) - z_mat(1:end-1);
                conv = sum(abs(real(diff(end-10:end))));
                if conv > 10
                    root = 0;
                else
                    root = 1;
                end
            else
                root = 1;
            end
            
            out = [z root];
            
        end
        
        function test_roots(obj,w,b)
            obj.omega = w;
            obj.num = obj.num;
            obj.eigenfunctionW(obj,w,b);
        end
        
        function value = vg(obj,w,b)
            
            dw      =   obj.eigenfunctionW(obj,w,b);
            dk      =   obj.eigenfunctionK(obj,w,b);
            
            value  =   -dk(2)/(obj.siC*dw(2));
        end
        
        function value = vgk(obj,w,b)
            
            dw      =   obj.eigenfunctionW(obj,w,b);
            dk      =   obj.eigenfunctionK(obj,w,b);
            
            tmp     =   -dw(2)/(dk(2));
            value   =   1/(obj.siC*real(tmp));
        end
        
        function value = vgk2(obj,w,b)
            
            dw      =   obj.eigenfunctionW(obj,w,b);
            dk      =   obj.eigenfunctionK(obj,w,b);
            
            tmp     =   -dw(2)/(dk(2));
            value   =   1/(obj.siC*tmp);
        end
        
        function value = ve(obj,w,b,xst,xfi,res)
            
            data = obj.energy_density(w,b,xst,xfi,res);
            
            vex = trapz(data(:,1), data(:,2))/ trapz(data(:,1), data(:,end));
            vey = trapz(data(:,1), data(:,3))/ trapz(data(:,1), data(:,end));
            vez = trapz(data(:,1), data(:,4))/ trapz(data(:,1), data(:,end));
            
            value = [vex/obj.siC vey/obj.siC vez/obj.siC];
            
        end
        
        function value = energy_density(obj,w,b,xst,xfi,res)
            
            data = obj.field_plotter(w,b,xst,xfi,res,0);
            field = data{1};
            y = flipud(field(:,1));
            fz = flipud(field(:,2));
            fy = flipud(field(:,3));
            fx = flipud(field(:,4));
            
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                Ex = zeros(length(fz),1);
                Ey = Ex;
                Ez = fz;
                
                Hx = fx;
                Hy = fy;
                Hz = Ex;
            else
                Hx = zeros(length(fz),1);
                Hy = Hx;
                Hz = fz;
                
                Ex = fx;
                Ey = fy;
                Ez = Hx;
            end
            
            for j = 1:length(data{2})
                if j == 1
                    eps = data{2}{j};
                    deps = data{3}{j};
                    mu = data{4}{j};
                    dmu = data{5}{j};
                else
                    eps = [eps data{2}{j}];
                    deps = [deps data{3}{j}];
                    mu = [mu data{4}{j}];
                    dmu = [dmu data{5}{j}];
                end
            end
            
            eps = fliplr(eps);
            deps = fliplr(deps);
            mu = fliplr(mu);
            dmu = fliplr(dmu);
            
            Ex2 = abs(Ex.^2);
            Ey2 = abs(Ey.^2);
            Ez2 = abs(Ez.^2);
            E2 = Ex2 + Ey2 + Ez2;
            
            Hx2 = abs(Hx.^2);
            Hy2 = abs(Hy.^2);
            Hz2 = abs(Hz.^2);
            H2 = Hx2 + Hy2 + Hz2;
            
            Sx = 0.5 * real(Ey.*conj(Hz) - Ez.*conj(Hy));
            Sy = 0.5 * real(Ez.*conj(Hx) - Ex.*conj(Hz));
            Sz = 0.5 * real(Ex.*conj(Hy) - Ey.*conj(Hx));
            
            EnEx = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ex2.' ;
            EnEy = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ey2.' ;
            EnEz = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ez2.' ;
            EnE = 0.25*obj.eps0*real(deps*obj.omega + eps).*E2.' ;
            
            EnHx = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hx2.';
            EnHy = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hy2.';
            EnHz = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hz2.';
            EnH = 0.25*obj.mu0*real(dmu*obj.omega + mu).*H2.';
            
            EnTot = 0.25*(obj.eps0*real(deps*obj.omega + eps).*E2.' + obj.mu0*(dmu*obj.omega + mu).*H2.');
            
            value = [y Sx Sy Sz EnEx.' EnEy.' EnEz.' EnE.' EnHx.' EnHy.' EnHz.' EnH.' EnTot.'];
        end
        
        function value = energy_density_RTA(obj,w,b,xst,xfi,res)
            
            data = obj.field_plotter_RTA(w,b,xst,xfi,res,0);
            field = data{1};
            y = flipud(field(:,1));
            fz = flipud(field(:,2));
            fy = flipud(field(:,3));
            fx = flipud(field(:,4));
            
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                Ex = zeros(length(fz),1);
                Ey = Ex;
                Ez = fz;
                
                Hx = fx;
                Hy = fy;
                Hz = Ex;
            else
                Hx = zeros(length(fz),1);
                Hy = Hx;
                Hz = fz;
                
                Ex = fx;
                Ey = fy;
                Ez = Hx;
            end
            
            for j = 1:length(data{2})
                if j == 1
                    eps = data{2}{j};
                    deps = data{3}{j};
                    mu = data{4}{j};
                    dmu = data{5}{j};
                else
                    eps = [eps data{2}{j}];
                    deps = [deps data{3}{j}];
                    mu = [mu data{4}{j}];
                    dmu = [dmu data{5}{j}];
                end
            end
            
            eps = fliplr(eps);
            deps = fliplr(deps);
            mu = fliplr(mu);
            dmu = fliplr(dmu);
            
            
            Ex2 = abs(Ex.^2);
            Ey2 = abs(Ey.^2);
            Ez2 = abs(Ez.^2);
            E2 = Ex2 + Ey2 + Ez2;
            
            
            
            Hx2 = abs(Hx.^2);
            Hy2 = abs(Hy.^2);
            Hz2 = abs(Hz.^2);
            H2 = Hx2 + Hy2 + Hz2;
            
            Sx = 0.5 * real(Ey.*conj(Hz) - Ez.*conj(Hy));
            Sy = 0.5 * real(Ez.*conj(Hx) - Ex.*conj(Hz));
            Sz = 0.5 * real(Ex.*conj(Hy) - Ey.*conj(Hx));
            

            EnEx = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ex2.' ;
            EnEy = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ey2.' ;
            EnEz = 0.25*obj.eps0*real(deps*obj.omega + eps).*Ez2.' ;
            EnE = 0.25*obj.eps0*real(deps*obj.omega + eps).*E2.' ;
            

            
            EnHx = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hx2.';
            EnHy = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hy2.';
            EnHz = 0.25*obj.mu0*real(dmu*obj.omega + mu).*Hz2.';
            EnH = 0.25*obj.mu0*real(dmu*obj.omega + mu).*H2.';
            
            EnTot = 0.25*(obj.eps0*real(deps*obj.omega + eps).*E2.' + obj.mu0*(dmu*obj.omega + mu).*H2.');
            
            
            value = [y Sx Sy Sz EnEx.' EnEy.' EnEz.' EnE.' EnHx.' EnHy.' EnHz.' EnH.' EnTot.'];
        end
        
        function value = field_plotter(obj,w,b,xst,xfi,res,pl)
            obj.omega = w;
            obj.num = obj.num;
            coe = obj.coefficients(obj,w,b);
            int_d = sum(obj.d);
            full_d = int_d+xst+xfi;
            inc = -full_d/res;
            n2_val      =   obj.n.^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val(2:end-1))-(b^2));
            
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val = obj.mu;
                mdiv = -obj.mu;
                fy = 'E_z'; %to match FDTD
                fx = 'H_y';
                fz = 'H_x';
                mul = obj.mu0;
            else
                m_val = obj.eps;
                mdiv = obj.eps;
                fy = 'H_z';
                fx = 'E_y';
                fz = 'E_x';
                mul = obj.eps0;
            end
            
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            
            gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
            gammaN      =   sgn(2)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
            
            xsup = xst:inc:0;
            
            epssup{1} = ones(1,length(xsup))*obj.eps(1);
            depssup{1} = ones(1,length(xsup))*obj.deps(1);
            musup{1} = ones(1,length(xsup))*obj.mu(1);
            dmusup{1} = ones(1,length(xsup))*obj.dmu(1);
            
            phisupy = exp(-gamma1*xsup);
            phisupx = b/(w*mdiv(1))*phisupy;
            phisupz = 1i*gamma1/(w*mdiv(1))*phisupy;
            
            numlay = length(obj.d);
            
            xmat = xsup.';
            phiymat = phisupy.';
            phixmat = phisupx.';
            phizmat = phisupz.';
            
            
            for loop = 1:numlay
                x = 0:inc:-obj.d(loop);
                epssup{1+loop} = ones(1,length(x))*obj.eps(1+loop);
                depssup{1+loop} = ones(1,length(x))*obj.deps(1+loop);
                musup{1+loop} = ones(1,length(x))*obj.mu(1+loop);
                dmusup{1+loop} = ones(1,length(x))*obj.dmu(1+loop);
                m = m_val(loop+1);
                phi1 = coe(1,loop);
                phi2 = coe(2,loop);
                ktm = kappa_val(loop);
                philayery = phi1*cos(ktm*x) + phi2*m/ktm*sin(ktm*x);
                philayerx = b/(w*mdiv(loop+1))*philayery;
                philayerz = 1i*ktm/(w*mdiv(loop+1))*(phi1*sin(ktm*x) - phi2*m/ktm*cos(ktm*x));
                phiymat = [phiymat; philayery.'];
                phixmat = [phixmat; philayerx.'];
                phizmat = [phizmat; philayerz.'];
                if loop == 1
                    xmat = [xmat; x.'];
                else
                    xmat = [xmat;x.'-sum(obj.d(1:loop-1))];
                end
            end
            
            
            xsub = 0:inc:-xfi;
            epssup{2+loop} = ones(1,length(xsub))*obj.eps(end);
            depssup{2+loop} = ones(1,length(xsub))*obj.deps(end);
            musup{2+loop} = ones(1,length(xsub))*obj.mu(end);
            dmusup{2+loop} = ones(1,length(xsub))*obj.dmu(end);
            phisuby = coe(1,end)*exp(gammaN*xsub);
            phisubx = b/(w*mdiv(end))*phisuby;
            phisubz = -1i*gammaN/(w*mdiv(end))*phisuby;
            
            xmat = [xmat;xsub.'-sum(obj.d)];
            phiymat = [phiymat; phisuby.'];
            phixmat = [phixmat; phisubx.']./mul;
            phizmat = [phizmat; phisubz.']./mul;
            
            if pl
                figure
                set(gcf,'color','w')
                subplot(2,2,1)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,real(phiymat),'linewidth',3);
                hold on
                %plot(xmat,imag(phiymat),'r--','linewidth',2);
                plot(xmat,zeros(1,length(xmat)),'k--');
                ax = axis;
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fy),'fontsize',20,'fontname','arial')
                
                subplot(2,2,2)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,real(phixmat),'linewidth',3);
                hold on
                %plot(xmat,imag(phixmat),'r--','linewidth',2);
                ax = axis;
                plot(xmat,zeros(1,length(xmat)),'k--');
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fx),'fontsize',20,'fontname','arial')
                
                subplot(2,2,3)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,real(phizmat),'linewidth',3);
                hold on
                %plot(xmat,imag(phizmat),'r--','linewidth',2);
                plot(xmat,zeros(1,length(xmat)),'k--');
                ax = axis;
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fz),'fontsize',20,'fontname','arial')
            end
            value = {[xmat phiymat phixmat phizmat],epssup,depssup,musup,dmusup};
        end
        
        function value = field_plotter_RTA(obj,w,b,xst,xfi,res,pl)
            
            obj.omega = w;
            obj.num = obj.num;
            tmp = obj.RTA_coefficients(obj,w,b);
            incident = tmp{1};
            trans = tmp{2};
            reflec = tmp{3};
            coe = tmp{4};
            int_d = sum(obj.d);
            full_d = int_d;
            inc = -full_d/res;
            n2_val      =   obj.n.^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val)-(b^2));
            
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val = obj.mu;
                mdiv = -obj.mu;
                fy = 'E_z'; %to match FDTD
                fx = 'H_y';
                fz = 'H_x';
                mul = obj.mu0;
            else
                m_val = obj.eps;
                mdiv = obj.eps;
                fy = 'H_z';
                fx = 'E_y';
                fz = 'E_x';
                mul = obj.eps0;
            end
            
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
                        
            xsup = xst:inc:0;
            
            epssup{1} = ones(1,length(xsup))*obj.eps(1);
            depssup{1} = ones(1,length(xsup))*obj.deps(1);
            musup{1} = ones(1,length(xsup))*obj.mu(1);
            dmusup{1} = ones(1,length(xsup))*obj.dmu(1);
            
            phisupy = incident*exp(1i*kappa_val(1)*xsup);
            phisupx = b/(w*mdiv(1))*phisupy;
            phisupz = kappa_val(1)/(w*mdiv(1))*phisupy;
            
            refsupy = reflec*exp(-1i*kappa_val(1)*xsup);
            refsupx = b/(w*mdiv(1))*refsupy;
            refsupz = kappa_val(1)/(w*mdiv(1))*refsupy;
                        
            numlay = length(obj.d);
            xmat = xsup.';
            phiymat = phisupy.';
            phixmat = phisupx.';
            phizmat = phisupz.';
            
            
            

            for loop = 1:numlay
                x = 0:inc:-obj.d(loop);
                epssup{loop+1} = ones(1,length(x))*obj.eps(1+loop);
                depssup{loop+1} = ones(1,length(x))*obj.deps(1+loop);
                musup{loop+1} = ones(1,length(x))*obj.mu(1+loop);
                dmusup{loop+1} = ones(1,length(x))*obj.dmu(1+loop);
                m = m_val(loop+1);
                phi1 = coe(1,loop);
                phi2 = coe(2,loop);
                ktm = kappa_val(loop+1);
                philayery = phi1*cos(ktm*x) + phi2*m/ktm*sin(ktm*x);
                philayerx = b/(w*mdiv(loop+1))*philayery;
                philayerz = 1i*ktm/(w*mdiv(loop+1))*(phi1*sin(ktm*x) - phi2*m/ktm*cos(ktm*x));
                phiymat = [phiymat; philayery.'];
                phixmat = [phixmat; philayerx.'];
                phizmat = [phizmat; philayerz.'];
                if loop == 1
                    xmat = [xmat; x.'];
                else
                    xmat = [xmat;x.'-sum(obj.d(1:loop-1))];
                end
            end
            
            
            
            xsub = 0:inc:-xfi;
            epssup{2+loop} = ones(1,length(xsub))*obj.eps(end);
            depssup{2+loop} = ones(1,length(xsub))*obj.deps(end);
            musup{2+loop} = ones(1,length(xsub))*obj.mu(end);
            dmusup{2+loop} = ones(1,length(xsub))*obj.dmu(end);
            
            phisuby = trans*exp(-1i*kappa_val(end)*xsub);
            phisubx = b/(w*mdiv(end))*phisuby;
            phisubz = -kappa_val(end)/(w*mdiv(end))*phisuby;
            
            xmat = [xmat;xsub.'-sum(obj.d)];
            phiymat = [phiymat; phisuby.'];
            phixmat = [phixmat; phisubx.']./mul;
            phizmat = [phizmat; phisubz.']./mul;
            refsupx = refsupx./mul;
            refsupz = refsupz./mul;
            
            if polac == 1

                E2 = abs(phiymat.^2);
                E2r = abs(refsupy.^2);
            else

                E2 = abs(phixmat.^2) + abs(phizmat.^2);
                E2r = abs(refsupx.^2) + abs(refsupz.^2);
            end
            
            if pl
                figure
                set(gcf,'color','w')
                subplot(2,2,1)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,abs(phiymat.^2),'linewidth',3);
                hold on
                plot(xsup,abs(refsupy.^2),'r','linewidth',3)
                %plot(xmat,imag(phiymat),'r--','linewidth',2);
                plot(xmat,zeros(1,length(xmat)),'k--');
                ax = axis;
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fy),'fontsize',20,'fontname','arial')
                
                subplot(2,2,2)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,abs(phixmat.^2),'linewidth',3);
                hold on
                plot(xsup,abs(refsupx.^2),'r','linewidth',3)
                
                %plot(xmat,imag(phixmat),'r--','linewidth',2);
                ax = axis;
                plot(xmat,zeros(1,length(xmat)),'k--');
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fx),'fontsize',20,'fontname','arial')
                
                subplot(2,2,3)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,abs(phizmat.^2),'linewidth',3);
                hold on
                plot(xsup,abs(refsupz.^2),'r','linewidth',3)
                
                %plot(xmat,imag(phizmat),'r--','linewidth',2);
                plot(xmat,zeros(1,length(xmat)),'k--');
                ax = axis;
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel(sprintf('%s (.arb)',fz),'fontsize',20,'fontname','arial')
                
                
                subplot(2,2,4)
                set(gca,'fontsize',16,'fontname','arial')
                plot(xmat,E2,'linewidth',3);
                hold on
                plot(xsup,E2r,'r','linewidth',3)
                
                %plot(xmat,imag(phizmat),'r--','linewidth',2);
                plot(xmat,zeros(1,length(xmat)),'k--');
                ax = axis;
                plot([0 0],[ax(3) ax(4)],'k')
                for loop = 1:length(obj.d)
                    inter = -sum(obj.d(1:loop));
                    plot([1 1]*inter,[ax(3) ax(4)],'k')
                end
                xlim([min(xmat) max(xmat)])
                ylim(ax(3:4))
                xlabel('y (m)','fontsize',20,'fontname','arial')
                ylabel('E2 (.arb)','fontsize',20,'fontname','arial')
            end
            value = {[xmat phiymat phixmat phizmat E2],epssup,depssup,musup,dmusup};
        end
        
    end
    
    methods (Static)
        % Materials ------------------------------------------------
        function value = Air(~,fid,i)
            fprintf(fid,'eps(%g) = 1;\n',i);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'deps(%g) = 0;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            value = [1 1 0 0];
        end
        
        function value = Si(~,fid,i)
            fprintf(fid,'eps(%g) = 11.68;\n',i);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'deps(%g) = 0;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            value = [11.68 1 0 0];
        end
        
        function value = Dielectric(eps,mu,fid,i)
            fprintf(fid,'eps(%g) = %g;\n',i,eps);
            fprintf(fid,'mu(%g) = %g;\n',i,mu);
            fprintf(fid,'deps(%g) = 0;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            value = [eps mu 0 0];
        end
        
        function value = Metal(eps_inf,wp,gamma,w,fid,i)
            tmp = Structure.Drude(eps_inf,wp,gamma,w,fid,i);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = NRI(eps_inf,wp,gamma,w,fid,i)
            tmp = Structure.DDrude(eps_inf,wp,gamma,w,fid,i);
            value = tmp;
        end
        
        function value = SiO2(fid,i)
            fprintf(fid,'eps(%g) = 3.46;\n',i);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'deps(%g) = 0;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            value = [3.46 1 0 0];
        end
        
        function value = ITO(w,fid,i)
            tmp = Structure.Drude(4,3.13e15,1.07e14,w,fid,i);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = ITO_4k(w,fid,i)
            tmp = Structure.Drude(4,3.13e15,1.07e12,w,fid,i);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = Andreas_drude(w,fid,i)
            tmp = Structure.Drude(1,1.39e16,3.2e13,w,fid,i);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = ITO_lossless(w,fid,i)
            tmp = Structure.Drude(4,3.13e15,0,w,fid,i);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = NRI_Hughes(w,fid,i)
            tmp = Structure.Drude(1,2*pi*490e12,2*pi*2e12,w,fid,i);
            m = Structure.Lorentz(1,2*pi*189.4e12,2*pi*165.4e12,2*pi*2e12,w,fid,i);
            value = [tmp(1) m(1) tmp(2) m(2)];
        end
        
        function value = NRI_Hughes_lossless(w,fid,i)
            e = Structure.Drude(1,2*pi*490e12,0,w,fid,i);
            m = Structure.Lorentz(1,2*pi*189.4e12,2*pi*165.4e12,0,w,fid,i);
            value = [e(1) m(1) e(2) m(2)];
        end
        
        function value = Drude(eps_inf,wp,gamma,w,fid,i)
            fprintf(fid,'eps(%g) = %g - (%g^2)./((w.^2)+complex(0,1).*%g.*w);\n',i,...
                eps_inf,wp,gamma);
            fprintf(fid,'deps(%g) = (2*w + (complex(0,%g)))*((%g^2)/((((w^2)+complex(0,1)*%g*w))^2));\n',i,...
                gamma,wp,gamma);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            drudeps = eps_inf - (wp^2)./((w.^2)+complex(0,1).*gamma.*w);
            depsdw = (2*w+(complex(0,gamma)))*((wp^2)/((((w^2)+complex(0,1)*gamma*w))^2));
            value = [drudeps depsdw];
        end
        
        function value = DDrude(eps_inf,wp,gamma,w,fid,i)
            fprintf(fid,'eps(%g) = %g - (%g^2)./((w.^2)+complex(0,1).*%g.*w);\n',i,...
                eps_inf,wp,gamma);
            fprintf(fid,'deps(%g) = (2*w + (complex(0,%g)))*((%g^2)/((((w^2)+complex(0,1)*%g*w))^2));\n',i,...
                gamma,wp,gamma);
            fprintf(fid,'mu(%g) = %g - (%g^2)./((w.^2)+complex(0,1).*%g.*w);\n',i,...
                eps_inf,wp,gamma);
            fprintf(fid,'dmu(%g) = (2*w + (complex(0,%g)))*((%g^2)/((((w^2)+complex(0,1)*%g*w))^2));\n',i,...
                gamma,wp,gamma);
            drudeps = eps_inf - (wp^2)./((w.^2)+complex(0,1).*gamma.*w);
            depsdw = (2*w+(complex(0,gamma)))*((wp^2)/((((w^2)+complex(0,1)*gamma*w))^2));
            value = [drudeps drudeps depsdw depsdw];
        end
        
        function value = Lorentz(eps_inf,mul,wp,gamma,w,fid,i)
            fprintf(fid,'eps(%g) = %g - %g*(%g^2)./((%g^2)-(w.^2)-2*complex(0,1).*%g.*w);\n',i,...
                eps_inf,mul,wp,wp,gamma);
            fprintf(fid,'deps(%g) = -2*%g*%g^2*(w + 1i*%g)/(((%g^2)-2*complex(0,1)*%g*w-(w^2))^2);\n',i,...
                mul,wp,gamma,wp,gamma);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            val = eps_inf - mul*(wp^2)./(wp^2-w.^2-complex(0,1).*gamma.*w);
            dval = -2*mul*wp^2*(1i*gamma + w)./((wp^2 - 1i*2*gamma*w - w.^2).^2);
            value = [val dval];
        end
        
        function value = four_lvl(obj,eps_inf,lam_a,lam_e,d_a,d_e,g_a,g_e,dens,inv,fid,i)
            
            w_e = 2*pi*299792458/lam_e;
            w_a = 2*pi*299792458/lam_a;
            s_e = obj.calc_sigma(w_e,d_e);
            s_a = obj.calc_sigma(w_a,d_a);
                   
            epssi = 8.8541878e-12;
            
            fprintf(fid,'n_e = (%g^2-w^2 + 1i*2*%g*w);\n',w_e,g_e);
            fprintf(fid,'d_e = ((%g^2-w^2)^2 + 4*%g^2*w^2);\n',w_e,g_e);
            fprintf(fid,'chi_e = (%g*%g*%g/%g) * n_e/d_e;\n',s_e,inv,dens,epssi);
            fprintf(fid,'n_a = (%g^2-w^2 + 1i*2*%g*w);\n',w_a,g_a);
            fprintf(fid,'d_a = ((%g^2-w^2)^2 + 4*%g^2*w^2);\n',w_a,g_a);
            fprintf(fid,'chi_a = (%g*%g*%g/%g) * n_a/d_a;\n',s_a,(1-inv),dens,epssi);
            fprintf(fid,'dn_e = -2*w + 1i*2*%g;\n',g_e);
            fprintf(fid,'dd_e = 2*w*(2*w^2 + 4*%g^2 - 2*%g^2);\n',g_e,w_e);
            fprintf(fid,'dchi_e = (%g*%g*%g/%g) * (dn_e*d_e - n_e*dd_e)/(d_e^2);\n',s_e,inv,dens,epssi);
            fprintf(fid,'dn_a = -2*w + 1i*2*%g;\n',g_a);
            fprintf(fid,'dd_a = 2*w*(2*w^2 + 4*%g^2 - 2*%g^2);\n',g_a,w_a);
            fprintf(fid,'dchi_a = (%g*%g*%g/%g) * (dn_a*d_a - n_a*dd_a)/(d_a^2);\n',s_a,(1-inv),dens,epssi);
            fprintf(fid,'eps(%g) = %g + chi_a - chi_e;\n',i,eps_inf);
            fprintf(fid,'deps(%g) = dchi_a - dchi_e;\n',i);
            fprintf(fid,'mu(%g) = 1;\n',i);
            fprintf(fid,'dmu(%g) = 0;\n',i);
            w = obj.omega;
            n_e = (w_e^2-w^2 + 1i*2*g_e*w);
            d_e = ((w_e^2-w^2)^2 + 4*g_e^2*w^2);
            chi_e = (s_e*inv*dens/epssi) * n_e/d_e;
            n_a = (w_a^2-w^2 + 1i*2*g_a*w);
            d_a = ((w_a^2-w^2)^2 + 4*g_a^2*w^2);
            chi_a = (s_a*(1-inv)*dens/epssi) * n_a/d_a;
            
            dn_e = -2*w + 1i*2*g_e;
            dd_e = 2*w*(2*w^2 + 4*g_e^2 - 2*w_e^2);
            dchi_e = (s_e*inv*dens/epssi) * (dn_e*d_e - n_e*dd_e)/(d_e^2);
            
            dn_a = -2*w + 1i*2*g_a;
            dd_a = 2*w*(2*w^2 + 4*g_a^2 - 2*w_a^2);
            dchi_a = (s_a*(1-inv)*dens/epssi)*(dn_a*d_a - n_a*dd_a)/(d_a^2);
            
            eps = eps_inf + chi_e - chi_a;
            deps = dchi_e - dchi_a;
            value = [eps deps 1 0];
            
            
            
        end
        
        function value = Silver_seb(w)
            eps_inf     =   1.17152;
            wd      =   1.39604e16;
            gammad  =   12.6126;
            dl1     =   2.23994;
            wl1     =   8.25718e15;
            gammal1 =   1.95614e14;
            dl2     =   0.222651;
            wl2     =   3.05707e15;
            gammal2 =   8.52675e14;
            
            i       =   complex(0,1);
            eps_s   =   eps_inf - (wd^2)/((w^2)-i*gammad*w);
            eps_s   =   eps_s - (dl1*(wl1^2))/((w^2)-(i*2*gammal1*w)-(wl1^2));
            eps_s   =   eps_s - (dl2*(wl2^2))/((w^2)-(i*2*gammal2*w)-(wl2^2));
            
            d_eps_s =   (wd^2)*(2*w-(i*gammad))/(((w^2)-i*gammad*w)^2);
            d_eps_s =   d_eps_s + (2*dl1*(wl1^2)*(w-i*gammal1))/(((w^2)-i*2*gammal1*w-(wl1^2))^2);
            d_eps_s =   d_eps_s + (2*dl2*(wl2^2)*(w-i*gammal2))/(((w^2)-i*2*gammal2*w-(wl2^2))^2);
            
            value     =   [eps_s 1 d_eps_s 0];
        end
        
        function value = Silver(w)
            tmp = Structure.Drude(4.03,1.39e16,3.14e13,w);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = Gold(w)
            tmp = Structure.Drude(9.0685,1.3544e16,1.1536e14,w);
            value = [tmp(1),1,tmp(2),0];
        end
        
        function value = Silver_lossless(w)
            tmp = Structure.Drude(1,1.39e16,0,w);
            value = [tmp(1),1,tmp(2),0];
        end
        % -----------------------------------------------------------
        
        % Dispersion calculations -----------------------------------
        
        function sol = eigenfunctionK(obj,w,b)
            
            n2_val      =   obj.n.^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val(2:end-1))-(b^2));
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val     =   obj.mu;
            else
                m_val     =   obj.eps;
            end
            
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            
            d_val       =   obj.d;
            gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
            gammaN      =   sgn(2)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
            
            % Initialise phi and psi
            phi1    =   1;
            phi2    =   -gamma1/m_val(1);
            phi     =   [phi1;phi2];
            
            psi1    =   0;
            psi2    =   -1/(2*m_val(1));
            psi     =   [psi1;psi2];
            
            N       =   length(m_val);
            M       =   zeros(2);
            
            for layer = 2:N-1
                
                dtm     =   d_val(layer-1);
                m       =   m_val(layer);
                kappa   =   kappa_val(layer-1);
                
                M(1,:)  =   [cos(kappa*dtm) (-(m/kappa)*sin(kappa*dtm))];
                M(2,:)  =   [(kappa/m)*sin(kappa*dtm) (cos(kappa*dtm))];
                
                DM0      =  (dtm/(2*kappa))*sin(kappa*dtm);
                DM1      =  ((-m/(2*(kappa^3)))*sin(kappa*dtm))+(m*dtm/(2*(kappa^2)))*cos(kappa*dtm);
                DM2      =  ((-1/(2*kappa*m))*sin(kappa*dtm))-((dtm/(2*m))*cos(kappa*dtm));
                DM3     =   (dtm/(2*kappa))*sin(kappa*dtm);
                
                DM      =   [DM0 DM1; DM2 DM3];
                
                psi     =   M*psi + gamma1*DM*phi;
                phi     =   M*phi;
                
            end
            
            H   =   phi(2)-(gammaN/m_val(end))*phi(1);
            DH  =   2*b*( (psi(2)/gamma1) - (phi(1)/(2*gammaN*m_val(end))) - (gammaN/(gamma1*m_val(end)))*psi(1));
            sol = [H DH];
        end
        
        function sol = eigenfunctionW(obj,w,b)
            
            n2_val      =   obj.n.^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val(2:end-1))-(b^2));
            dkappa_val = (obj.n(2:end-1)*w).*(obj.dn(2:end-1)*w+...
                obj.n(2:end-1))./(kappa_val*obj.siC^2);
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val = obj.mu;
                dm_val = obj.dmu;
            else
                m_val = obj.eps;
                dm_val = obj.deps;
            end
            
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            
            d_val       =   obj.d;
            gamma1      =   sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
            dgamma1     =   -w*obj.n(1)*(w*obj.dn(1)+obj.n(1))/(gamma1*(obj.siC^2));
            gammaN      =   sqrt((b^2)-((w/obj.siC)^2)*n2_val(end));
            dgammaN     =   -w*obj.n(end)*(w*obj.dn(end)+obj.n(end))/(gammaN*(obj.siC^2));
            gamma1      =   sgn(1)*gamma1;
            dgamma1     =   sgn(1)*dgamma1;
            gammaN      =   sgn(2)*gammaN;
            dgammaN     =   sgn(2)*dgammaN;
            
            
            % Initialise phi and psi
            phi1    =   1;
            phi2    =   -gamma1/m_val(1);
            phi     =   [phi1;phi2];
            
            psi1    =   0;
            psi2    =   ((gamma1*dm_val(1))-(dgamma1*m_val(1)))/(m_val(1)^2);
            psi     =   [psi1;psi2];
            
            N       =   length(m_val);
            M = zeros(2);
            
            for layer = 2:N-1
                
                dtm       =   d_val(layer-1);
                m       =   m_val(layer);
                dm      =   dm_val(layer);
                kappa   =   kappa_val(layer-1);
                d_k     =   dkappa_val(layer-1);
                
                M(1,:)       =   [cos(kappa*dtm) (-(m/kappa)*sin(kappa*dtm))];
                M(2,:)       =   [(kappa/m)*sin(kappa*dtm) (cos(kappa*dtm))];
                
                DM0     =   -dtm*d_k*sin(dtm*kappa);
                DM1     =   (((m*d_k)-(dm*kappa))/(kappa^2))*sin(kappa*dtm);
                DM1     =   DM1 - ((m*dtm*d_k/kappa)*cos(kappa*dtm));
                DM2     =   (((m*d_k)-(dm*kappa))/(m^2))*sin(kappa*dtm);
                DM2     =   DM2 + ((dtm*kappa*d_k/m)*cos(kappa*dtm));
                
                DM      =   [DM0 DM1; DM2 DM0];
                
                psi     =   M*psi + DM*phi;
                phi     =   M*phi;
                
            end
            
            H   =   phi(2)-(gammaN/m_val(end))*phi(1);
            
            DH  =   (psi(2)) + (((gammaN*dm_val(end))-(dgammaN*m_val(end)))/(m_val(end)^2))*phi(1)-((gammaN*psi(1))/(m_val(end)));
            
            sol = [H DH];
        end
        
        function sol = coefficients(obj,w,b)
            n2_val      =   obj.n.^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val(2:end-1))-(b^2));
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val = obj.mu;
            else
                m_val = obj.eps;
            end
            
            switch obj.bnd
                case 0
                    sgn = [1 1];
                case 1
                    sgn = [-1 1];
                case 2
                    sgn = [1 -1];
                case 3
                    sgn = [-1 -1];
                otherwise
                    fprintf('Please choose bound = 0-3');
                    sgn = [1 1];
            end
            
            d_val       =   obj.d;
            gamma1      =   sgn(1)*sqrt((b^2)-((w/obj.siC)^2)*n2_val(1));
            
            % Initialise phi and psi
            phi1    =   1;
            phi2    =   -gamma1/m_val(1);
            phi     =   [phi1;phi2];
            
            N       =   length(m_val);
            M = zeros(2);
            phimat = zeros(2,N-1);
            phimat(:,1) = phi;
            for layer = 2:N-1
                dtm     =   d_val(layer-1);
                m       =   m_val(layer);
                kappa   =   kappa_val(layer-1);
                
                M(1,:)       =   [cos(kappa*dtm) (-(m/kappa)*sin(kappa*dtm))];
                M(2,:)       =   [(kappa/m)*sin(kappa*dtm) (cos(kappa*dtm))];
                
                phi     =   M*phi;
                phimat(:,layer) = phi;
            end
            
            sol = phimat;
        end
        
        function sol = RTA(obj,w,b_in)
            
            obj.omega = w;
            obj.beta = b_in;
            obj.num = obj.num;
            c1      =   [0;1];
            c_par   =   c1;
            c_per   =   c1;
            k0 = w/obj.siC;
            nsinc   =   b_in/k0;
            n_struc = obj.n;
            z_struc = obj.impz;
            d_struc = obj.d;
            
            % Flip structure around
            n_struc = fliplr(n_struc);
            z_struc = fliplr(z_struc);
            d_struc = fliplr(d_struc);
            d_struc = [d_struc 0];
            
            
            % Calculate amplitude coefficients
            
            for loop2 = 1:length(n_struc)-1
                
                n1      =   n_struc(loop2);
                n2      =   n_struc(loop2+1);
                stheta1 =   nsinc/n1;
                stheta2 =   nsinc/n2;
                ctheta1 =   sqrt(1-stheta1^2);
                ctheta2 =   sqrt(1-stheta2^2);
                z1      =   z_struc(loop2);
                z2      =   z_struc(loop2+1);
                d2      =   d_struc(loop2);
                
                a       =   ctheta1/ctheta2;
                b       =   z2/z1;
                s       =   n2*k0*ctheta2;
                e_pos   =   exp(complex(0,-1)*s*d2);
                e_neg   =   exp(complex(0,1)*s*d2);
                
                M_par0  =   ((a+b)*e_neg);
                M_par1  =   ((a-b)*e_neg);
                M_par2  =   ((a-b)*e_pos);
                M_par3  =   ((a+b)*e_pos);
                
                M_par   =   0.5.*[M_par0 M_par1;M_par2 M_par3];
                
                M_per0  =   (1+(a*b))*e_neg;
                M_per1  =   (1-(a*b))*e_neg;
                M_per2  =   (1-(a*b))*e_pos;
                M_per3  =   (1+(a*b))*e_pos;
                
                M_per   =   0.5.*[M_per0 M_per1;M_per2 M_per3];
                
                c_par   =   M_par*c_par;
                c_per   =   M_per*c_per;
                
            end
            
            % Calculate Power coefficients
            
            sthetac =   nsinc/n_struc(1);
            sthetas =   nsinc/n_struc(end);
            
            cthetac =   sqrt(1-sthetac^2);
            cthetas =   sqrt(1-sthetas^2);
            
            zs      =   z_struc(end);
            zc      =   z_struc(1);
            inc_par =   (c_par(2));
            tau_par =   1/inc_par;
            ref_par =   c_par(1)/inc_par;
            R_par   =   abs(ref_par)^2;
            
            
            T_par   =   (zs/zc)*(cthetac/cthetas)*abs(tau_par)^2;
            
            inc_per =   c_per(2);
            tau_per =   1/inc_per;
            ref_per =   c_per(1)/inc_per;
            R_per   =   abs(ref_per)^2;
            
            T_per   =   (zs/zc)*(cthetac/cthetas)*abs(tau_per)^2;
            A_par       =   1-(R_par+T_par);
            A_per       =   1-(R_per+T_per);
            sol = [R_par T_par A_par R_per T_per A_per];
            
        end
        
        function sol = RTA_New(obj,w,b)
            
            obj.omega = w;
            obj.beta = b;
            obj.num = obj.num;
            n2_val      =   fliplr(obj.n).^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val)-(b^2));
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val     =   fliplr(obj.mu);
            else
                m_val     =   fliplr(obj.eps);
            end
            
            d_val       =   fliplr(obj.d);
            
            
            % Initialise phi and psi
            phi1    =   1;
            phi2    =   1i*kappa_val(1)/m_val(1);
            phi     =   [phi1;phi2];
            
            N       =   length(m_val);
            M       =   zeros(2);
            
            for layer = 2:N-1
                
                dtm     =   d_val(layer-1);
                m       =   m_val(layer);
                
                kappa   =   kappa_val(layer);
                
                M(1,:)  =   [cos(kappa*dtm) (-(m/kappa)*sin(kappa*dtm))];
                M(2,:)  =   [(kappa/m)*sin(kappa*dtm) (cos(kappa*dtm))];
                
                phi     =   M*phi;
                
            end
            
            B1 = (0.5*(phi1 - 1i*phi2*m_val(1)/kappa_val(1)));
            An = (0.5*(phi(1) + 1i*phi(2)*m_val(end)/kappa_val(end)));
            Bn = (0.5*(phi(1) - 1i*phi(2)*m_val(end)/kappa_val(end)));
            
            Pi = conj(kappa_val(end))/(2*w*conj(m_val(end))) * conj(Bn)*Bn;
            Pr = conj(kappa_val(end))/(2*w*conj(m_val(end))) * conj(An)*An;
            Pt = conj(kappa_val(1))/(2*w*conj(m_val(1))) * conj(B1) * B1;
            
            T = Pt/Pi;
            R = Pr/Pi;
            A = 1-(R+T);
            
            sol = [R T A];
            
        end
        
        function sol = RTA_coefficients(obj,w,b)
            
            obj.omega = w;
            obj.beta = b;
            obj.num = obj.num;
            n2_val      =   fliplr(obj.n).^2;
            kappa_val   =   sqrt((((w/obj.siC)^2)*n2_val)-(b^2));
            polac = isequal(obj.pola,'TE');
            
            if polac == 1
                m_val     =   fliplr(obj.mu);
            else
                m_val     =   fliplr(obj.eps);
            end
            
            d_val       =   fliplr(obj.d);
            
            
            % Initialise phi and psi
            phi1    =   1;
            phi2    =   -1i*kappa_val(1)/m_val(1);
            phi     =   [phi1;phi2];
            
            N       =   length(m_val);
            M       =   zeros(2);
            phimat(:,1) = phi;
            for layer = 2:N-1
                
                dtm     =   -d_val(layer-1);
                m       =   m_val(layer);
                
                kappa   =   kappa_val(layer);
                
                M(1,:)  =   [cos(kappa*dtm) (-(m/kappa)*sin(kappa*dtm))];
                M(2,:)  =   [(kappa/m)*sin(kappa*dtm) (cos(kappa*dtm))];
                
                phi     =   M*phi;
                
                phimat(:,layer) = phi;
                
            end
            
            B1 = (0.5*(phi1 + 1i*phi2*m_val(1)/kappa_val(1)));
            Bn = (0.5*(phi(1) + 1i*phi(2)*m_val(end)/kappa_val(end)));
            An = (0.5*(phi(1) - 1i*phi(2)*m_val(end)/kappa_val(end)));
            sol = {Bn,B1,An,fliplr(phimat)};
            
        end
        
        % -----------------------------------------------------------
        
        % Miscelanious ----------------------------------------------
        function value = roundd(num,idp)
            numin   =   num;
            
            for lp = 1:length(num)
                num     =   real(numin(lp));
                if num == 0
                    reout = 0;
                else
                    s       =   sign(num);
                    num     =   s.*num;
                    mul     =   floor(log10(num));
                    num     =   num.*10.^-mul;
                    mult    =   10^idp;
                    tmp     =   floor((num*mult)+0.5)/mult;
                    reout(lp,1)     =   s.*tmp.*10.^mul;
                end
                
                num     =   imag(numin(lp));
                if num == 0
                    imout = 0;
                else
                    s       =   sign(num);
                    num     =   s.*num;
                    mul     =   floor(log10(num));
                    num     =   num.*10.^-mul;
                    mult    =   10^idp;
                    tmp     =   floor((num*mult)+0.5)/mult;
                    imout(lp,1)     =   s.*tmp.*10.^mul;
                end
            end
            
            value = complex(reout,imout);
        end
        
        function value = pow3(input)
            p = floor(log(input)/log(1000));
            ls = {'y','z','a','f','p','n','u','m',...
                '','k','M','G','T','P','E','Z','Y'};
            value = {1000^p, p, ls{p+9}};
        end
        
        function value = closest(x,x0)
            ind = find(abs(x-x0) == min(abs(x-x0)));
            value = ind(1);
        end
        
        function output = calc_sigma(w,d)
            eV     =   1.6021765e-19;      % Elementary charge (C)
            hbar   =   1.0545717e-34;      % Reduced Planck's constant (Js)
            sigma    =   (2*w * eV^2 * d^2) / ...
                (hbar);
            
            output  =   sigma;
            
        end
        % -----------------------------------------------------------
    end
    
end %Classdef