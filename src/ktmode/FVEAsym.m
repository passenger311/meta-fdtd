

lamda=1.55e-6; 
 hx=dx*1e-6; 
 hy=dy*1e-6; 
% 
% W=2e-6; %4e-6; %3e-6;
% ws=1e-6; %0.99e-6;
% ena=0.5e-6; %1e-6; ;
% dio=0.52e-6; %0.30e-6;
% tria=2e-6; %1.35e-6; 
% %tese=0.1e-6; %0.15e-6;
% %pente=0.245e-6; %0.1e-6;
% %eksi=0.1e-6;
% %epta=0.05e-6;
% keh=1e-6;
% kew=0.5e-6; %0.5e-6;
% 
% %w=W/2+kew;
% w=W/2+kew;
% 
% z1=round(ws/2/hx);
% z2=round(W/2/hx);
%PREVIOUS BY KOSMAS kenowr=129; %round(w/hx);

% enar=round(ena/hy);
% dior=round(dio/hy)+enar;
% triar=round(tria/hy)+dior;
%teser=round(tese/hy)+triar;
%penter=round(pente/hy)+teser;
%eksir=round(eksi/hy)+penter;
%eptar=round(epta/hy)+eksir;

%PREVIOUS BY KOSMAS kenoh2r=158; %round(keh/hy)+triar;
%kenoh2r=round(keh/hy)+eptar;

%DANIELE DEFINES kenoh2r,kenowr
[kenoh2r,kenowr] = size(N);

tot=kenowr*kenoh2r;
tot2=2*tot;

%n=zeros(1,tot);

Lamminus1xx=zeros(1,tot);
Lamxx=zeros(1,tot);
Lamplus1xx=zeros(1,tot);
Subdiag1xx=zeros(1,tot-kenowr);
Subdiag2xx=zeros(1,tot-1);
Diagxx=zeros(1,tot);
Subdiag3xx=zeros(1,tot-1);
Subdiag4xx=zeros(1,tot-kenowr);

Subdiag1xy=zeros(1,tot-kenowr-1);
Subdiag2xy=zeros(1,tot-kenowr+1);
Subdiag3xy=zeros(1,tot-kenowr+1);
Subdiag4xy=zeros(1,tot-kenowr-1);
Lamxy1=zeros(1,tot);
Lamxy2=zeros(1,tot);
Lamxy3=zeros(1,tot);
Lamxy4=zeros(1,tot);


%- - - - - - - - - - - - - - - - - 


Lamminus1yy=zeros(1,tot);
Lamyy=zeros(1,tot);
Lamplus1yy=zeros(1,tot);
Subdiag1yy=zeros(1,tot-kenowr);
Subdiag2yy=zeros(1,tot-1);
Diagyy=zeros(1,tot);
Subdiag3yy=zeros(1,tot-1);
Subdiag4yy=zeros(1,tot-kenowr);

Subdiag1yx=zeros(1,tot-kenowr-1);
Subdiag2yx=zeros(1,tot-kenowr+1);
Subdiag3yx=zeros(1,tot-kenowr+1);
Subdiag4yx=zeros(1,tot-kenowr-1);
Lamyx1=zeros(1,tot);
Lamyx2=zeros(1,tot);
Lamyx3=zeros(1,tot);
Lamyx4=zeros(1,tot);

B=zeros(tot2,13);

for j=1:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Lamminus1xx(ii)=2*(n(ii-1)^2)/((n(ii-1)^2)+(n(ii)^2));
    end
end

for j=1:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamplus1xx(ii)=2*(n(ii+1)^2)/((n(ii)^2)+(n(ii+1))^2);
    end
end

for j=1:kenoh2r
    for i=2:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamxx(ii)=2+((n(ii)^2)-(n(ii+1)^2))/((n(ii)^2)+(n(ii+1)^2))+((n(ii)^2)-(n(ii-1)^2))/((n(ii)^2)+(n(ii-1)^2));
    end
end

for j=1:kenoh2r
    ii=(j-1)*kenowr+1;
    Lamxx(ii)=2+((n(ii)^2)-(n(ii+1)^2))/((n(ii)^2)+(n(ii+1)^2));
end

for j=1:kenoh2r
    ii=(j-1)*kenowr+kenowr;
    Lamxx(ii)=2+((n(ii)^2)-(n(ii-1)^2))/((n(ii)^2)+(n(ii-1)^2));
end

Subdiag1xx(1,1:1:tot-kenowr)=1/(hy*hy);

for j=1:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag2xx(1,ii-1)=Lamminus1xx(ii)/(hx*hx);
    end
end

Diagxx(1,1:1:tot)=((2*pi*n(1,1:1:tot)/lamda).^2)-(Lamxx(1,1:1:tot)./(hx*hx))-(2/(hy*hy));

for j=1:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag3xx(1,ii)=Lamplus1xx(ii)/(hx*hx);
    end
end

% for j=1:kenoh2r
%     ii=(j-1)*kenowr+1;
%     Subdiag3xx(1,ii)=2*Lamplus1xx(ii)/(hx*hx);
% end

Subdiag4xx(1,1:1:tot-kenowr)=1/(hy*hy);

for j=2:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Lamxy1(ii)=((n(ii-kenowr-1)^2)/(n(ii-1)^2))-1;
    end
end

for j=2:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamxy2(ii)=-(((n(ii-kenowr+1)^2)/(n(ii+1)^2))-1);
    end
end

for j=1:kenoh2r-1
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Lamxy3(ii)=-(((n(ii+kenowr-1)^2)/(n(ii-1)^2))-1);
    end
end

for j=1:kenoh2r-1
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamxy4(ii)=((n(ii+kenowr+1)^2)/(n(ii+1)^2))-1;
    end
end

for j=2:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag1xy(1,ii-kenowr-1)=Lamxy1(ii)/(4*hx*hy);
    end
end

for j=2:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag2xy(1,ii-kenowr+1)=Lamxy2(ii)/(4*hx*hy);
    end
end

% for j=2:kenoh2r
%     ii=(j-1)*kenowr+1;
%     Subdiag2xy(1,ii-kenowr+1)=2*Lamxy2(ii)/(4*hx*hy);
% end

for j=1:kenoh2r-1
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag3xy(1,ii)=Lamxy3(ii)/(4*hx*hy);
    end
end

for j=1:kenoh2r-1
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag4xy(1,ii)=Lamxy4(ii)/(4*hx*hy);
    end
end

% for j=1:kenoh2r-1
%     ii=(j-1)*kenowr+1;
%     Subdiag4xy(1,ii)=2*Lamxy4(ii)/(4*hx*hy);
% end

for j=2:kenoh2r
    for i=1:kenowr
        ii=(j-1)*kenowr+i;
        Lamminus1yy(ii)=2*(n(ii-kenowr)^2)/((n(ii-kenowr)^2)+(n(ii)^2));
    end
end

for j=1:kenoh2r-1
    for i=1:kenowr
        ii=(j-1)*kenowr+i;
        Lamplus1yy(ii)=2*(n(ii+kenowr)^2)/((n(ii)^2)+(n(ii+kenowr))^2);
    end
end

for j=2:kenoh2r-1
    for i=1:kenowr
        ii=(j-1)*kenowr+i;
        Lamyy(ii)=2+((n(ii)^2)-(n(ii+kenowr)^2))/((n(ii)^2)+(n(ii+kenowr)^2))+((n(ii)^2)-(n(ii-kenowr)^2))/((n(ii)^2)+(n(ii-kenowr)^2));
    end
end

for i=1:kenowr
    Lamyy(i)=2+((n(i)^2)-(n(i+kenowr)^2))/((n(i)^2)+(n(i+kenowr)^2));
end

for i=1:kenowr
    ii=(kenoh2r-1)*kenowr+i;
    Lamyy(ii)=2+((n(ii)^2)-(n(ii-kenowr)^2))/((n(ii)^2)+(n(ii-kenowr)^2));
end

Subdiag1yy(1,1:1:tot-kenowr)=Lamminus1yy(1,kenowr+1:1:tot)./(hy*hy);

for j=1:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag2yy(1,ii-1)=1/(hx*hx);
    end
end

Diagyy(1,1:1:tot)=((2*pi*n(1,1:1:tot)/lamda).^2)-(Lamyy(1,1:1:tot)./(hy*hy))-(2/(hx*hx));

% for j=1:kenoh2r
%     ii=(j-1)*kenowr+1;
%     Diagyy(1,ii)=0;
% end

for j=1:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag3yy(1,ii)=1/(hx*hx);
    end
end

% for j=1:kenoh2r
%     ii=(j-1)*kenowr+1;
%     Subdiag3yy(1,ii)=0;
% end

Subdiag4yy(1,1:1:tot-kenowr)=Lamplus1yy(1,1:1:tot-kenowr)./(hy*hy);

for j=2:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Lamyx1(ii)=((n(ii-kenowr-1)^2)/(n(ii-kenowr)^2))-1;
    end
end

for j=2:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamyx2(ii)=-(((n(ii-kenowr+1)^2)/(n(ii-kenowr)^2))-1);
    end
end

for j=1:kenoh2r-1
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Lamyx3(ii)=-(((n(ii+kenowr-1)^2)/(n(ii+kenowr)^2))-1);
    end
end

for j=1:kenoh2r-1
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Lamyx4(ii)=((n(ii+kenowr+1)^2)/(n(ii+kenowr)^2))-1;
    end
end

for j=2:kenoh2r
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag1yx(1,ii-kenowr-1)=Lamyx1(ii)/(4*hx*hy);
    end
end

for j=2:kenoh2r
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag2yx(1,ii-kenowr+1)=Lamyx2(ii)/(4*hx*hy);
    end
end

% for j=2:kenoh2r
%     ii=(j-1)*kenowr+1;
%     Subdiag2yx(1,ii-kenowr+1)=0;
% end

for j=1:kenoh2r-1
    for i=2:kenowr
        ii=(j-1)*kenowr+i;
        Subdiag3yx(1,ii)=Lamyx3(ii)/(4*hx*hy);
    end
end

for j=1:kenoh2r-1
    for i=1:kenowr-1
        ii=(j-1)*kenowr+i;
        Subdiag4yx(1,ii)=Lamyx4(ii)/(4*hx*hy);
    end
end

% for j=1:kenoh2r-1
%     ii=(j-1)*kenowr+1;
%     Subdiag4yx(1,ii)=0;
% end

for i=1:tot-kenowr-1
    B(i,1)=Subdiag1yx(1,i);
    B(i+kenowr+1,4)=Subdiag4yx(1,i);
    B(i+tot,10)=Subdiag1xy(1,i);
    B(i+tot+kenowr+1,13)=Subdiag4xy(1,i);
end

for i=1:tot-kenowr+1
    B(i,2)=Subdiag2yx(1,i);
    B(i+kenowr-1,3)=Subdiag3yx(1,i);
    B(i+tot,11)=Subdiag2xy(1,i);
    B(i+tot+kenowr-1,12)=Subdiag3xy(1,i);
end

for i=1:tot-kenowr
    B(i,5)=Subdiag1xx(1,i);
    B(i+kenowr,9)=Subdiag4xx(1,i);
end

for i=1:tot-1
    B(i,6)=Subdiag2xx(1,i);
    B(i+1,8)=Subdiag3xx(1,i);
end

for i=1:tot
    B(i,7)=Diagxx(1,i);
end

for i=1:tot-kenowr
    B(i+tot,5)=Subdiag1yy(1,i);
    B(i+tot+kenowr,9)=Subdiag4yy(1,i);
end

for i=1:tot-1
    B(i+tot,6)=Subdiag2yy(1,i);
    B(i+tot+1,8)=Subdiag3yy(1,i);
end

for i=1:tot
    B(i+tot,7)=Diagyy(1,i);
end

d=[-(tot+kenowr+1);
   -(tot+kenowr-1);
   -(tot-kenowr+1);
   -(tot-kenowr-1);
    -kenowr;
    -1;
     0;
     1;
    kenowr;
    tot-kenowr-1;
    tot-kenowr+1;
    tot+kenowr-1;
    tot+kenowr+1];
    

S=spdiags(B,d,tot2,tot2);



