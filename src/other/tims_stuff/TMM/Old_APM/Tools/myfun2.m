
% -------------------------------------------------------------------------
% |                                                                       |
% |                           Function 'myfun2'                           |
% |                                                                       |
% |  This function is used in the function z_value and is the function    |
% |  that is integrated in the second part of the Arguement Pricinple     |
% |  Method equal to z*H'(z)/H(z).  The numerical integrator used feeds   |
% |  the function an array of input values along the line of integration. |
% |  As the function eigenfunction cant handle arrays of z values this    |
% |  function contains a for loop that feeds the eigenfunction one value  |
% |  at a time and then builds an array of results (f) to send back to    |
% |  the numerical integrator.                                            |
% -------------------------------------------------------------------------


function f = myfun2(z)

global R                            % Coefficient used in mapping
global S                            % Coefficient used in mapping
global s_par                        % Free space wavenumber
global p_par                        % Layer thickness array
global pola

Num = length(z);                    % Number of points to evaluate
f = [];                             % Initialize results array            

k0          =   p_par(4);           % Free space wavevector (1/m)
eps_val     =   s_par(1,:);
mu_val      =   s_par(2,:);
d_val       =   s_par(5,3:end);
if pola == 1
    m_val = s_par(2,:);
elseif pola == 0;
    m_val = s_par(1,:);
end
for loop = 1:Num
    % Call eigenfunction to calculate H and H' at z
    sol = eigenfunction(eps_val,d_val,m_val,mu_val,R,S,z(loop),k0);
    f1 = z(loop)*(sol(2)/sol(1));           % Calculate function z*H'/H
    f = [f f1];                             % Add result to array

end


end