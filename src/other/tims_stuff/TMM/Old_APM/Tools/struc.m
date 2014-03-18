
% -------------------------------------------------------------------------
% |                                                                       |
% |                            Function Struc                             |
% |                                                                       |
% | In this function you define the layered structure you want to solve.  |
% | The structure is always of the form shown below.  When solving for    |
% | reflection and transmission coefficients light is always incident     |
% | from the cover layer.                                                 |
% |                                                                       |               
% |                                                                       |
% |                      |                        |                       |                         
% |          Cover layer | .... Layers (1:N) .... | Substrate layer       |
% |             (epsc)   |          eps(1:N)      |     (epss)            |
% |             (muc)    |           mu(1:N)      |     (mus)             |
% |                      |            d(1:N)      |                       |
% |                                                                       |
% -------------------------------------------------------------------------

function par = struc(struc_name,p_par)

fh  =   str2func(struc_name);
par =   fh(p_par);