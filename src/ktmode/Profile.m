
%tic
format long
ko=2*pi/lamda;
%sh=(nco*ko)^2;
sh=1.8e14; %(20e6)^2; %8e13; %(nwa*ko)^2;
kenowrf=2*kenowr-1;
L1=zeros(kenoh2r,kenowr);
L2=zeros(kenoh2r,kenowr);

ndeik1=zeros(kenoh2r,kenowr);

% x=(-kenowr:1:kenowr-2)*hx*1e6;
% y=(-enar:1:kenoh2r-1-enar)*hy*1e6;
% [X,Y]=meshgrid(x,y);

[V,idio]=eigs(S,1,sh);
%[V,d]=eigs(S,1,'LR');

for j=1:kenoh2r
    for i=1:kenowr
        ii=(j-1)*kenowr+i;
        L1(j,i)=abs(V(ii));
        ndeik1(j,i)=abs(n(1,ii));
    end
end

for j=kenoh2r+1:2*kenoh2r
    for i=1:kenowr
        ii=(j-1)*kenowr+i;
        L2(j-kenoh2r,i)=abs(V(ii));
    end
end

a=max(abs(L1));
b1=max(a);
[imax,jmax]=find(abs(L1)==b1);
if L1(imax,jmax)<0
    b1=-b1;
%else
%    b=1;
end

%b2=max(a);
%[imax,jmax]=find(abs(L1)==b2);
%if L1(imax,jmax)<0%a=max(abs(L2));

%    b2=-b2;
%else
%    b=1;
%end

for j=1:kenowr
    for i=1:kenoh2r
        L1(i,j)=L1(i,j)/(b1);
    end
end

for j=1:kenowr
    for i=1:kenoh2r
        L2(i,j)=L2(i,j)/(b1);
    end
end

%figure

subplot(3,1,1)
surf(abs(L1));
axis tight;
shading interp;
colorbar;
view(0,-90);
% 
subplot(3,1,2)
surf(abs(L2));
axis tight;
shading interp;
colorbar;
view(0,-90);
% %axis([-kenowr*hx*1e6 (kenowr-1)*hx*1e6 -enar*hy*1e6 (kenoh2r-1-enar)*hy*1e6]);
% 
subplot(3,1,3)
surf(ndeik1);
axis tight;
shading interp;
colorbar;
view(0,-90);
% %axis([1 kenowrf 1 kenoh2r]);
% %axis([-kenowr*hx*1e6 (kenowr-1)*hx*1e6 -enar*hy*1e6 (kenoh2r-1-enar)*hy*1e6]);

bita=sqrt(idio);
neff=bita/ko;
%toc

