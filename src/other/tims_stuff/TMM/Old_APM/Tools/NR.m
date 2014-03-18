
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function 'NR'                               |
% |                                                                       |
% |  This function is performs the Newton Raphson method to find roots    |
% |  with greater accuracy and is also used in the root tracking part of  |
% |  the program.                                                         |
% -------------------------------------------------------------------------

function z = NR(eps_val,d_val,m_val,mu_val,R,S,z,k0,nrs)

% Perform Newton Raphson method for nrs number of times
for loop = 1:nrs

    sol = eigenfunction(eps_val,d_val,m_val,mu_val,R,S,z,k0);

    H = sol(1);
    DH = sol(2);

    z = z - (H/DH);    % Newton Raphson equation

end

end
