%function output = drude(w)
clear all
close all
c = 299792458;
loop=0;
for lambda = 300:1:2200
    loop = loop+1;
    w = (2*pi*c)/(lambda*1e-9);
wp          =   3.13e15;
gamma       =   1.07e14;
eps_inf     =   4;
eps         =   eps_inf-((wp^2)/((w^2)+complex(0,gamma*w)));
eps_mat(loop) = eps;
lambda_mat(loop) = lambda;
end
figure(1)

plot(lambda_mat,real(eps_mat),'r','LineWidth',2);
hold on
plot([300 3500],[0 0],'k');
hold on
plot([500 3500],[-15.7797 -15.7797],'k');
hold on
plot([1718 1718],[-20 10],'k');
hold on
plot([2981 2981],[-20 10],'k');
axis([300 2200 -10 4])
figure(2)
plot(lambda_mat,imag(eps_mat));
hold on
plot([1550 1550],[0 7]);
hold on
plot([2250 2250],[0 7]);
axis([300 2200 0 4])
%output = eps;
prec = 16;
% Save root results
FileName    =   'AZO_eps.txt';
FileID      =   fopen(FileName, 'w');

sym = '%';


    a = sprintf('lambda\tRe(eps)\tIm(eps)\n');
    fprintf(FileID,a);
    a1 = sprintf('%s.%ge\t%s+.%ge\t%s+.%ge\n',...
        sym,prec,sym,prec,sym,prec);
    n = size(lambda_mat,2);
    for loop = 1:n
        l       =   lambda_mat(loop);
        re       =   real(eps_mat(loop));
        im     =   imag(eps_mat(loop));

        
       
        fprintf(FileID,a1,l,re,im);
        
    end
    
    
% Close Text Files
fclose(FileID);