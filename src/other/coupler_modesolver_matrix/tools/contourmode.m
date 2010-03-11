function [cs,h] = contourmode(x,y,mode,maxm,dbsep,dbmax,t);

% Produces a contour plot (in dB) of the mode of an optical
% waveguide.  
% 
% USAGE:
% 
% contourmode(x,y,mode,maxm,dbsep,dbmax,t);
% [cs,h] = contourmode(x,y,mode,maxm,dbsep,dbmax,t);
% 
% INPUT:
% 
% x,y - vectors describing horizontal and vertical grid points
% mode - the mode or field component to be plotted
% maxm - the field amplitude that should correspond to 0 dB.
% dbsep - separation between successive contours (in dB)
% dbmax - the smallest contour to plot (in dB)
% t - title of plot
% 
% OUTPUT:
% 
% cs - contour matrix
% h - handle to contourgroup object

x = real(x);
y = real(y);
mode = abs(transpose(mode));
v = (0:-dbsep:-dbmax)';
[cs,h] = contour(x,y,20*log10(mode/maxm),v);
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
xlabel('x'); 
ylabel('y'); 
title(t);
v = axis;
set(gca,'PlotBoxAspectRatio',[v(2)-v(1) v(4)-v(3) 1]);
