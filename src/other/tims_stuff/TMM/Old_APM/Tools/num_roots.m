
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function 'num_roots'                        |
% |                                                                       |
% |  This function performs the contour integral of the Argument          |
% |  Principle method. It performs a numerical integration around the     |
% |  contour defined by the points re0,re1,im0,im1 in a counter-clockwise |
% |  direction:                                                           |
% |                                   3                                   |
% |                   (re0,im1)-------<<------(re1,im1)                   |
% |                       |                  .    |                       |
% |   Im ^                |    .                  |                       |
% |      |                |           .           |                       |
% |      |                v      .                ^                       |
% |      |              4 v                 .     ^ 2                     |
% |       ------>         |        (roots)        |                       |
% |             Re        |            .          |                       |
% |                       |    .                  |                       |
% |                       |                    .  |                       |
% |                   (re0,im0)------->>------(re1,im0)                   |
% |                                   1                                   |
% |                                                                       |
% |  using Matlabs inbuilt numerical integrator quadgk although other     |
% |  integrators such as quad or quadl can be used. The function attempts |
% |  to return the number of roots inside the contour with a minimum      |
% |  error but when the integration is difficult, particularly when there |
% |  is a large number of roots, the error might be large and the step    |
% |  size is reduced to attempt to gain greater accuracy however this     |
% |  will decrease performance.  The function myfun returns the value of  |
% |  (H'/H) used in the APM.                                              |
% -------------------------------------------------------------------------



function Num_R = num_roots(contour)

global tol

z(1)    =   complex(contour(1),contour(3));
z(2)    =   complex(contour(2),contour(3));
z(3)    =   complex(contour(2),contour(4));
z(4)    =   complex(contour(1),contour(4));
z(5)    =   z(1);

Q       =   zeros(1,4);


for n = 1:4
    Q(n)        =   quadgk(@myfun,z(n),z(n+1));
    Q_real      =   log(myfun3(z(n+1)))-log(myfun3(z(n)));
    Q_error     =   abs(real(Q(n))-real(Q_real));
    
    if Q_error > tol
        Q(n) = 0;
        a       =   z(n);
        b_fin   =   z(n+1);
        b       =   a+((b_fin-a)/2);
        
        
        while a ~= z(n+1)
            
            while Q_error > tol
                Q_par       =   quadgk(@myfun,a,b);
                Q_real      =   log(myfun3(b))-log(myfun3(a));
                Q_error     =   abs(real(Q_par)-real(Q_real));
                b_out       =   b;
                b           =   a+((b-a)/2);  
                
            end
            Q_error = tol*2;
            a = b_out;
            b = b_fin;
            Q(n) = Q(n)+Q_par;

            
        end
        

    end
    
end

num = (sum(Q))/(complex(0,2*pi));  % Total contour integral divided by i2pi

Num_R = round(real(num));             % Number of roots found

end