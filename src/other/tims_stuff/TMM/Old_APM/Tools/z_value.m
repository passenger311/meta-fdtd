
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function 'z_value'                          |
% |                                                                       |
% |  This function performs the second contour integral of the Arguement  |
% |  Principle Method z*H'(z)/H(z). It is used once a single root has     |
% |  been isolated and returns the position of the root.  Like the        |
% |  function num_roots it performs a numerical integration around the    |
% |  contour defined by the points re0,re1,im0,im1 in a counter-clockwise |
% |  direction:                                                           |
% |                                   3                                   |
% |                   (re0,im1)-------<<------(re1,im1)                   |
% |                       |                       |                       |
% |   Im ^                |                       |                       |
% |      |                |                       |                       |
% |      |                v     .                 ^                       |
% |      |              4 v                       ^ 2                     |
% |       ------>         |         (root)        |                       |
% |             Re        |                       |                       |
% |                       |                       |                       |
% |                       |                       |                       |
% |                   (re0,im0)------->>------(re1,im0)                   |
% |                                   1                                   |
% |                                                                       |
% |  using Matlabs inbuilt numerical integrator quadgk although other     |
% |  integrators such as quad or quadl can be used. The function attempts |
% |  to return the number of roots inside the contour with a minimum      |
% |  error but when the integration is difficult, when the root is close  |
% |  to the contour, the error might be large and the step size is        |
% |  reduced to attempt to gain greater accuracy however this will        |
% |  decrease performance.  The function myfun2 returns the value of      |
% |  z*(H/H') used in the APM.                                            |
% -------------------------------------------------------------------------



function z = z_value(re0,re1,im0,im1)

Q1 = quadgk(@myfun2,complex(re0,im0),complex(re1,im0)); % Int along path 1
Q2 = quadgk(@myfun2,complex(re1,im0),complex(re1,im1)); % Int along path 2
Q3 = quadgk(@myfun2,complex(re1,im1),complex(re0,im1)); % Int along path 3
Q4 = quadgk(@myfun2,complex(re0,im1),complex(re0,im0)); % Int along path 4

Q = (Q1+Q2+Q3+Q4)/(2*pi*i);  % Total contour integral divided by i2pi
z = Q;                       % is equal to the position (z) of the root   
end