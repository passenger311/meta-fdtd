function output = vg(w,b,pola,c,name,d_struc)

dw      =   dD_dw(w,b,pola,c,name,d_struc);
dk      =   dD_dk(w,b,pola,c,name,d_struc);

output  =   -real(dk(2))/(real((c*dw(2))));