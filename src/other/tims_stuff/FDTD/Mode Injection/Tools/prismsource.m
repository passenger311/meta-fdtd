% This program is designed to take the input file source-prism.txt and to
% construct a source file.

clear all
close all
%if exist('/user/phpgr/ph31ek/fdtd3d-par/trunk/Prism/source-prism2.txt') == 2
%ustop = 2
%else
%ustop = 1
%end

%for u = 1:ustop
%if u == 1 
workdir = cd;
a   =   ls('source_profile*');
a   =   textscan(a,'%s');
a   =   a{1};


Num_s = size(a,1);
for loop = 1:Num_s
    name = sprintf('source_profile%g.txt',loop);
file_input = fopen(name);
%elseif u == 2
%file_input = fopen('/home/atiraid/atitacpg/ph31ek/fdtd3d-par/trunk/Prism/source-prism2.txt');
%end 

A = textscan(file_input,'%s %n','headerlines',0);
        sourcexst     = A{2}(1,1);
        sourcexlength = A{2}(2,1);
        sourceyst     = A{2}(3,1);
        sourceyfi     = A{2}(4,1);
        sup_epsilon   = A{2}(5,1);
        sup_mu        = A{2}(6,1);
        sigma         = A{2}(7,1);



fclose(file_input);
        
 %NumberPoints = (sourcexfi-sourcexst)+1+(sourceyfi-sourceyst);
NumberPoints = sourcexlength+(sourceyfi-sourceyst); 

% E = zeros(1,NumberPoints);
% H = zeros(1,NumberPoints);
Eamp = sqrt(sup_mu/sup_epsilon);
Hamp = 1;


b = ceil(NumberPoints/2);
q = 1:NumberPoints;
E = Eamp*(exp((-1*(q-b).^2)/(2*sigma^2)));
H = Hamp*(exp((-1*(q-b).^2)/(2*sigma^2)));

% E = sqrt(sup_mu/sup_epsilon)*ones(1,NumberPoints);
% H = ones(1,NumberPoints);
%      
%      for q = 1:NumberPoints;
%          E(q) = E(q)*(exp((-1*(q-b)^2)/(2*sigma^2)));
%          H(q) = H(q)*(exp((-1*(q-b)^2)/(2*sigma^2)));
%      end

 
%if u == 1 
% file_output = fopen('/home/atiraid/atitacpg/ph31ek/fdtd3d-par/trunk/Prism/source-prism-data.txt','wt');
%elseif u == 2
% file_output = fopen('/home/atiraid/atitacpg/ph31ek/fdtd3d-par/trunk/Prism/source-prism-data2.txt','wt');
%end
%a = '/home/atiraid/atitacpg/ph31ek/fdtd3d-par/trunk/Prism/source-prism-data'
a = sprintf('%s/source_prism_data%g.in',workdir,loop)

L=sourcexlength;


file_output = fopen(a,'wt');

 fprintf(file_output,'(SET\n');
 fprintf(file_output,'%5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f\n',sourcexst, sourcexst+sourcexlength-1, 1, sourceyst, sourceyfi, 1, 0, 0, 1);
 
 for p = 1:NumberPoints;
     fprintf(file_output,'%12.6E %12.6E\n',E(p), H(p));
 end;
 
 fprintf(file_output,')SET');
fclose(file_output);

 q=(1:NumberPoints);
%  figure(1);plot(q,E,q,H,'--')
%  axis([q(1) q(end) 0 1])
pause(2)
end
quit;
