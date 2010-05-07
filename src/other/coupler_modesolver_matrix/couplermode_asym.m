function [neff]=couplermode(lambda,neffguess,nmodes,dx,dy,x,y,xc,yc,eps)


figure(2)
imagemode(xc,yc,log(real(eps))/log(max(max(abs(eps)))),'Epsilon structure');


fprintf (1,'solving for eigenmodes...'); t = cputime;
[Hx,Hy,neff] = wgmodes (lambda, neffguess, nmodes, dx, dy, eps, '0000');
fprintf (1,'done (cputime = %7.3f)\n', cputime-t);



if nmodes==1,
   fprintf (1,'post-processing...'); t = cputime;
   [Hz,Ex,Ey,Ez] = postprocess (lambda, neff, Hx, Hy, dx, dy, eps, '0000');
   fprintf (1,'done (cputime = %7.3f)\n', cputime-t);

   ii = 1;
   colormap(jet(256));
   hn = max(abs([Hy(:);Hx(:)]));
   %en = max(abs([Ey(:);Ex(:)]));
   en = hn/neff;
   %en = hn/ncore;
   figure(1);

   if 1 == 1
   
   subplot(231);
   imagemode(x,y,Hx/hn,sprintf('Hx (mode %d)',ii));
   hold on;
   contourmode(x,y,Hx/hn,1,3,45,sprintf('Hx (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(234);
   imagemode(x,y,Hy/hn,sprintf('Hy (mode %d)',ii));
   hold on;
   contourmode(x,y,Hy/hn,1,3,45,sprintf('Hy (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(233);
   imagemode(x,y,Hz/hn,sprintf('Hz (mode %d)',ii));
   hold on;
   contourmode(x,y,Hz/hn,1,3,45,sprintf('Hz (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);
   
   subplot(235);
   
   end;
   
   imagemode(xc,yc,Ex/en,sprintf('Ex (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ex/en,1,3,60,sprintf('Ex (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   if 1 == 1
   
   subplot(232);
   imagemode(xc,yc,Ey/en,sprintf('Ey (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ey/en,1,3,60,sprintf('Ey (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);

   subplot(236);
   imagemode(xc,yc,Ez/en,sprintf('Ez (mode %d)',ii));
   hold on;
   contourmode(xc,yc,Ez/en,1,3,60,sprintf('Ez (mode %d)',ii));
   hold off;
   v = xlim();
   %line(v,[h1,h1]);
   %line([0,w/2,w/2],[h1+h2,h1+h2,h1]);
   figure(2)
   imagemode(xc,yc,log(real(eps))/log(max(max(abs(eps)))),'Epsilon structure');

   end;
   
   end;
