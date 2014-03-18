%% Intialise Structure
close all
clear all
clc
test = Structure;
test.layers = {            % Define Layers
    'Air',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'SiO2',...
    'Si',...
    'Air'
    };
test.d = [113.4 269.1 113.4 269.1 113.4 269.1 113.4 269.1 200 269.1 113.4 269.1 113.4 269.1 113.4 269.1 113.4]*1e-9;
test.pola = 'TM';
test.bnd = 3;
ang = 76.444;

bin = 1.2e6; %2*pi/1.55e-6*cosd(ang);

%%
tic;
roots = test.find_roots_w(1,3e15,200,bin);
t = toc;

% Print Time taken and roots found
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));
msg = sprintf('Time taken to find %g roots: %gh %gm %.3gs\n',length(roots),floor(h),floor(m),s);
fprintf(1,msg);

for n = 1:length(roots)
    fprintf('\nRoot %g has omega = %g and loss = %g',n,real(roots(n)),...
        imag(roots(n))*2);
end

% Plot field profile
close all
for loop = 1:length(roots)
    
field = test.field_plotter(roots(loop),bin,500e-9,500e-9,2000);

end

% Ask for roots to track
a = input('\n\nPlease select root to track: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end
close all

data = test.geodispersionW(roots(a),bin,490e-9,510e-9,100,1);

svname = sprintf('Changing_ITO_d_76_degrees.txt');
dlmwrite(svname, data, 'delimiter', '\t', ...
         'precision', 13);
return

%% real dispersion

disp = test.realdispersion(0,2e7,100,0,3e15,100);
pcolor(disp{1},disp{2},log(abs(disp{3})))
shading interp
colormap(hot(1024))

%% reflection
reflecdata = test.reflec2D(1,3e15,200,1,2e7,100,'Beta');
x = reflecdata{1};
y = reflecdata{2};
figure(1)
pcolor(x,y,reflecdata{3})
shading interp
figure(2)
pcolor(x,y,real(reflecdata{4}))
shading interp
figure(3)
pcolor(x,y,real(reflecdata{5}))
shading interp
figure(4)
pcolor(x,y,real(reflecdata{6}))
shading interp
figure(5)
pcolor(x,y,real(reflecdata{7}))
shading interp
figure(6)
pcolor(x,y,real(reflecdata{7}))
shading interp
return
%% Find Complex w roots
tic;
roots = test.find_roots_w(1e15,1.5e15,200,bin);
t = toc;

% Print Time taken and roots found
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));
msg = sprintf('Time taken to find %g roots: %gh %gm %.3gs\n',length(roots),floor(h),floor(m),s);
fprintf(1,msg);

for n = 1:length(roots)
    fprintf('\nRoot %g has omega = %g and loss = %g',n,real(roots(n)),...
        imag(roots(n))*2);
end

% Plot field profile
close all
for loop = 1:length(roots)
    
field = test.field_plotter(roots(loop),bin,500e-9,500e-9,200);

end

% Ask for roots to track
a = input('\n\nPlease select root to track: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end
close all

% Find zgv points
% zgv = test.findvp(bin,roots(5),0)
% return

% Track roots
tic
if matlabpool('size') > 0
    matlabpool close force local;
end

NCores = feature('NumCores');
%NCores = 2;
matlabpool(NCores);

parfor (i = 1:length(a),NCores)
    fprintf(1,sprintf('\nTracking mode %g of %g\n',i,length(a)));
    mode{i} = test.dispersionW(roots(a(i)),bin,0,2e7,200);
    loss_mat(:,i) = mode{i}(:,5);
    k_mat(:,i) = mode{i}(:,3);
    vg_mat(:,i) = mode{i}(:,6);
    str{i} = sprintf('Mode %g',a(i));
end

t = toc;
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));



% for n = 1:length(a)
%     mode{n} = test.dispersionW(roots(a(n)),bin,0,5e7,1000);
%     loss_mat(:,n) = mode{n}(:,5);
%     k_mat(:,n) = mode{n}(:,3);
%     vg_mat(:,n) = mode{n}(:,6);
%     str{n} = sprintf('Mode %g',a(n));
% end
% toc

% Plot Roots
figure
set(gcf,'color','w')
plot(mode{1}(:,1),k_mat,'LineWidth',3)
hold on
plot(mode{1}(:,1),mode{1}(:,1)*test.siC,'k--')
plot(mode{1}(:,1),mode{1}(:,1)*test.siC/sqrt(11.68),'k--')
set(gca,'FontSize',28,'FontName','Arial')
xlabel('\beta, (m^{-1})','FontSize',30,'FontName','Arial')
ylabel('\omega, (rad s^{-1})','FontSize',30,'FontName','Arial')
ylim([0 2e15])
legend(str,'Location','Best')

figure
set(gcf,'color','w')
semilogy(mode{1}(:,1),abs(real(vg_mat)),'LineWidth',3)
set(gca,'FontSize',28,'FontName','Arial')
xlabel('\beta, (m^{-1})','FontSize',30,'FontName','Arial')
ylabel('v_g, (c)','FontSize',30,'FontName','Arial')
legend(str,'Location','Best')


figure
set(gcf,'color','w')
plot(mode{1}(:,1),-loss_mat,'LineWidth',3)
set(gca,'FontSize',28,'FontName','Arial')
xlabel('\beta, (m^{-1})','FontSize',30,'FontName','Arial')
ylabel('Loss, \alpha (s^{-1})','FontSize',30,'FontName','Arial')
legend(str,'Location','Best')



for loop = 1:length(mode)
tmp = mode{loop};

M = [tmp(:,1:5),real(tmp(:,6)),sign(tmp(:,7:end))];

svname = sprintf('Mode_coreplusplus_%g.txt',loop);
dlmwrite(svname, M, 'delimiter', '\t', ...
         'precision', 13);
end

msg = sprintf('\n\nTime taken to track roots: %gh %gm %.3gs\n',floor(h),floor(m),s);
fprintf(1,msg);
return

%% Find complex k roots
win = 1.21e15;
tic;
roots = test.find_roots_k(1e6,6e6,100,win);
t = toc;

% Print Time taken and roots found
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));
msg = sprintf('Time taken to find %g roots: %gh %gm %.3gs\n',length(roots),floor(h),floor(m),s);
fprintf(1,msg);

for n = 1:length(roots)
    fprintf('\nRoot %g has beta = %g and loss = %g',n,real(roots(n)),...
        imag(roots(n))*2);
end

% Plot field profile
close all
for loop = 1:length(roots)
    
field = test.field_plotter(win,roots(loop),500e-9,500e-9,2000);

end

% Ask for roots to track
a = input('\n\nPlease select root to track: ','s');
if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end

% Track roots
for n = 1:length(a)
    mode{n} = test.dispersionK(roots(a(n)),win,1.207e15,1.217e15,1000);
    str{n} = sprintf('Mode %g',a(n));
    loss_mat(:,n) = mode{n}(:,5);
    k_mat(:,n) = mode{n}(:,3);
    w_mat(:,n) = mode{n}(:,1);
    vg_mat(:,n) = mode{n}(:,6);
end

% Plot Roots
figure
for FLOOP = 1:2
    plot(k_mat(:,FLOOP),w_mat(:,FLOOP),'LineWidth',3)
    hold on
end
legend(str)
return
w = mode{1}(:,1);
k1 = k_mat(:,2);
difw = w(2:end)-w(1:end-1);
difk = k1(2:end)-k1(1:end-1);
vg = difw./difk;
vg = vg/299792458;
w2 = (w(2:end)+w(1:end-1))/2;
k2 = (k1(2:end)+k1(1:end-1))/2;
figure
plot(k2,real(vg))
legend(str)
t = toc;
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));
% 
% km = mode{1}(:,1);
% om = k_mat;
% 
% 
% ind = find(km>=(om/test.siC));
tmp = mode{1};

M = [tmp(:,1:5),real(tmp(:,6)),sign(tmp(:,7:end))];

svname = 'Complex_k_left_zoom.txt';
dlmwrite(svname, M, 'delimiter', '\t', ...
         'precision', 8);

msg = sprintf('\n\nTime taken to track roots: %gh %gm %.3gs\n',floor(h),floor(m),s);
fprintf(1,msg);




%%
%plot(refl(:,1),refl(:,3))
return
disp = test.realdispersion(1,2e7,200,1,3e15,200);
pcolor(disp{1},disp{2},log(abs(disp{3})))
shading interp
colormap(hot(1024))
t = toc;
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));
msg = sprintf('Time taken to find %g roots: %gh %gm %.3gs\n',length(roots),floor(h),floor(m),s);
fprintf(1,msg);

for n = 1:length(roots)
    fprintf('\nRoot %g has omega = %g and loss = %g',n,real(roots(n)),...
        imag(roots(n))*2);
end

a = input('\n\nPlease select root to track: ','s');
tic

if strcmpi(a(1),'a')
    a = 1:length(roots);
else
    a = str2num(a);
end

for n = 1:length(a)
    mode{n} = test.dispersionK(roots(a(n)),1.5e15,0,3e15,1000);
    k_mat(:,n) = mode{n}(:,3);
    str{n} = sprintf('Mode %g',a(n));
end


plot(k_mat,mode{1}(:,1))

legend(str)
t = toc;
h = t/60/60;
m = 60*(h-floor(h));
s = 60*(m-floor(m));

msg = sprintf('Time taken to track root: %gh %gm %.3gs\n',floor(h),floor(m),s);
fprintf(1,msg);
